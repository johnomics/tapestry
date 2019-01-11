import pysam

from collections import namedtuple
from statistics import mean

from Bio.SeqUtils import GC

from .misc import memoize

# Defining contig report at top level rather than in class so it works with multiprocessing
def contig_report(contig):
    tel_start, tel_end = contig.num_telomeres()
    return f"{contig.name}\t{len(contig)}\t{contig.GC()}\t{contig.mean_depth('reads')}\t{contig.mean_depth('contigs')}\t{contig.num_alignments()}\t{tel_start}\t{tel_end}"
            
DepthRecord = namedtuple('DepthRecord', 'start, end, depth')

class Contig:
    def __init__(self, rec, telomeres, outdir):
        self.name = rec.id
        self.rec = rec
        self.telomeres = telomeres
        self.outdir = outdir

    def __len__(self):
        return len(self.rec.seq)

    @memoize
    def GC(self):
       return f"{GC(self.rec.seq):.1f}"

    @memoize
    def depths(self, mapped):
        depths=[]
        if os.path.exists(f"{self.outdir}/{mapped}_assembly.regions.bed.gz"):
            for line in zgrep(f"^{self.name}", f"{self.outdir}/{mapped}_assembly.regions.bed.gz").split('\n')[:-1]: # Lose final empty line with :-1
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
