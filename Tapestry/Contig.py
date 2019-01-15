import os, pysam

from collections import namedtuple
from statistics import mean

from intervaltree import Interval, IntervalTree

from Bio.SeqUtils import GC

from .misc import memoize, grep

# Defining assembly reports at top level rather than in class so they work with multiprocessing
def contig_report(contig):
    tel_start, tel_end = contig.num_telomeres()
    report = f"{contig.name}"
    report += f"\t{len(contig)}"
    report += f"\t{contig.GC()}"
    report += f"\t{contig.mean_depth('reads')}"
    report += f"\t{contig.mean_depth('contigs')}"
    report += f"\t{contig.num_alignments()}"
    report += f"\t{tel_start}"
    report += f"\t{tel_end}"
    report += f"\t{contig.pc_unique()}"

    return report


def redundancy_report(contig):
    regions = ""
    for region in contig.region_depths():
        regions += f"{contig.name}\t{region.begin}\t{region.end}\t{region.end-region.begin}\t{region.data}\n"
    return regions.rstrip()


DepthRecord = namedtuple('DepthRecord', 'start, end, depth')

class PAF:
    def __init__(self, pafline):
        f = pafline.rstrip().split('\t')
        self.query_name      = f[0]
        self.query_length    = int(f[1])
        self.query_start     = int(f[2])
        self.query_end       = int(f[3])
        self.strand          = f[4]
        self.target_name     = f[5]
        self.target_length   = int(f[6])
        self.target_start    = int(f[7])
        self.target_end      = int(f[8])
        self.matches         = int(f[9])
        self.block_length    = int(f[10])
        self.mapping_quality = int(f[11])


class Contig:
    def __init__(self, rec, telomeres, outdir):
        self.name = rec.id
        self.rec = rec
        self.telomeres = telomeres
        self.outdir = outdir

    def __len__(self):
        return len(self.rec.seq)

    def __lt__(self, other):
        return len(self) < len(other)

    @memoize
    def GC(self):
       return f"{GC(self.rec.seq):.1f}"

    @memoize
    def depths(self, mapped):
        depths=[]
        if os.path.exists(f"{self.outdir}/{mapped}_assembly.regions.bed.gz"):
            for line in grep(f"^{self.name}", f"{self.outdir}/{mapped}_assembly.regions.bed.gz"): # Lose final empty line with :-1
                contigname, start, end, depth = line.split('\t')
                depths.append(DepthRecord(start=int(start), end=int(end), depth=float(depth))) 
        return depths

    @memoize
    def mean_depth(self, mapped):
        return f"{mean([d.depth for d in self.depths(mapped)]):.1f}" if self.depths(mapped) else 0

    @memoize
    def num_alignments(self):
        alignments=0
        if os.path.exists(f"{self.outdir}/reads_assembly.bam"):
            bam = pysam.AlignmentFile(f"{self.outdir}/reads_assembly.bam", 'rb')
            for aln in bam.fetch(self.name):
                alignments += 1
        return alignments

    @memoize
    def num_telomeres(self):
        start_matches = end_matches = 0
        if self.telomeres:
            for t in self.telomeres:
                for s in t, t.reverse_complement():
                    start_matches += len(list(s.instances.search(self.rec[:1000].seq)))
                    end_matches   += len(list(s.instances.search(self.rec[-1000:].seq)))
        return start_matches, end_matches

    @memoize
    def contig_alignments(self):
        alignments = IntervalTree()
        alignments[1:len(self)] = 1
        if os.path.exists(f"{self.outdir}/contigs_assembly.paf.gz"):
            for line in grep(f"{self.name}", f"{self.outdir}/contigs_assembly.paf.gz"):
                alignment = PAF(line)
                if alignment.query_name == self.name:
                    alignments[alignment.query_start:alignment.query_end] = 1
                if alignment.target_name == self.name:
                    alignments[alignment.target_start:alignment.target_end] = 1
        return alignments

    @memoize
    def region_depths(self):
        alignments = self.contig_alignments()
        regions = alignments.copy()
        regions.split_overlaps()
        region_depths = IntervalTree()
        for region in regions:
            region_depths[region.begin:region.end] = len(alignments[region.begin:region.end])
        
        return sorted(region_depths)

    @memoize
    def pc_unique(self):
        unique_bases = len(self)
        for region in self.region_depths():
            if region.data > 1:
                unique_bases -= region.end - region.begin # No need to -1 because end is beyond upper limit
        return f"{unique_bases/len(self) * 100:.0f}"