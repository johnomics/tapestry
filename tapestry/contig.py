import os, pysam

import numpy as np

from collections import namedtuple, defaultdict
from statistics import mean, median

from intervaltree import Interval, IntervalTree

from sklearn import mixture

from Bio.SeqUtils import GC

from .misc import grep, PAF

# Define process_contig at top level rather than in class so it works with multiprocessing
def process_contig(contig):
    contig.process()
    return contig

def get_ploidy(contig, median_depth=None):
    contig.ploidys = contig.get_ploidys(median_depth)
    return contig

DepthRecord = namedtuple('DepthRecord', 'start, end, depth')


class Contig:

    def __init__(self, rec, telomeres, outdir):
        self.name = rec.id
        self.rec = rec
        self.telomeres = telomeres
        self.outdir = outdir


    def report(self, assembly_gc):
        report = f"{self.name}"
        report += f"\t{len(self)}"
        report += f"\t{self.gc:.1f}"
        report += f"\t{self.median_read_depth:.1f}"
        report += f"\t{self.mean_read_depth:.1f}"
        report += f"\t{self.mean_contig_depth:.1f}"
        report += f"\t{self.tel_start}"
        report += f"\t{self.tel_end}"
        report += f"\t{self.mean_start_overhang}"
        report += f"\t{self.mean_end_overhang}"
        report += f"\t{self.unique_bases}"
        report += f"\t{self.unique_pc:.0f}"
        report += f"\t{self.category(assembly_gc)}"
        for p in sorted(self.ploidys):
            report += f"\t{p}:{self.ploidys[p]:.2f}"
    
        return report


    def __len__(self):
        return len(self.rec.seq)


    def __lt__(self, other):
        return len(self) < len(other)


    def category(self, assembly_gc):
        category=''
        if abs(self.gc - assembly_gc) < 2:
            category += 'N'
        else:
            category += '-'

        completeness = ''
        if self.tel_start > 0 and self.mean_start_overhang < 250:
            completeness += 'L'
        if self.tel_end   > 0 and self.mean_end_overhang   < 250:
            completeness += 'R'
        if completeness == 'LR':
            completeness = 'C'
        category += completeness if completeness else '-'

        max_ploidy = sorted(self.ploidys, key=lambda x:self.ploidys[x], reverse=True)[0] if self.ploidys else 0
        if max_ploidy > 4:
            max_ploidy = 'R'
        category += f"{max_ploidy}"

        return category


    def redundancy_report(self):
        regions = ""
        for region in self.region_depths:
            regions += f"{self.name}\t{region.begin}\t{region.end}\t{region.end-region.begin}\t{region.data}\n"
        return regions.rstrip()


    def process(self):
        self.gc = self.get_gc()
        self.contig_depths = self.depths('contigs')
        self.read_depths = self.depths('reads')
        self.mean_contig_depth = self.mean_depth(self.contig_depths)
        self.mean_read_depth = self.mean_depth(self.read_depths)
        self.median_read_depth = self.median_depth(self.read_depths)
        self.contig_alignments = self.get_contig_alignments()
        self.mean_start_overhang, self.mean_end_overhang = self.get_read_overhangs()
        self.region_depths = self.get_region_depths()
        self.unique_bases = self.get_unique_bases()
        self.unique_pc = self.get_unique_pc()
        self.tel_start, self.tel_end = self.num_telomeres()


    def get_gc(self):
        return GC(self.rec.seq)


    def depths(self, mapped):
        depths=[]
        if os.path.exists(f"{self.outdir}/{mapped}_assembly.regions.bed.gz"):
            for line in grep(f"^{self.name}", f"{self.outdir}/{mapped}_assembly.regions.bed.gz"): # Lose final empty line with :-1
                contigname, start, end, depth = line.split('\t')
                depths.append(DepthRecord(start=int(start), end=int(end), depth=float(depth))) 
        return depths


    def mean_depth(self, depths):
        return mean([d.depth for d in depths]) if depths else 0


    def median_depth(self, depths):
        depths = [d.depth for d in depths]
        return median(depths) if depths else 0


    def get_read_overhang(self, bam, start, end, overhang_function):
        num_reads = overhang_bases = 0
        for aln in bam.fetch(self.name, start, end):
            num_reads += 1
            overhang = overhang_function(aln)
            if overhang > 0:
                overhang_bases += overhang
        mean_overhang = int(overhang_bases / num_reads) if num_reads > 0 else None
        return mean_overhang


    def get_read_overhangs(self):
        mean_start_overhang = mean_end_overhang = 0
        if os.path.exists(f"{self.outdir}/reads_assembly.bam"):
            bam = pysam.AlignmentFile(f"{self.outdir}/reads_assembly.bam", 'rb')

            mean_start_overhang = self.get_read_overhang(
                    bam, 0, 1000,
                    lambda aln: aln.query_alignment_start-aln.reference_start)
            mean_end_overhang = self.get_read_overhang(
                    bam, len(self)-1001, len(self)-1,
                    lambda aln: (aln.query_length - aln.query_alignment_end) - (len(self) - aln.reference_end))

        return mean_start_overhang, mean_end_overhang


    def num_telomeres(self):
        start_matches = end_matches = 0
        if self.telomeres:
            for t in self.telomeres:
                for s in t, t.reverse_complement():
                    start_matches += len(list(s.instances.search(self.rec[:1000].seq)))
                    end_matches   += len(list(s.instances.search(self.rec[-1000:].seq)))
        return start_matches, end_matches


    def get_contig_alignments(self):
        alignments = IntervalTree()
        alignments[1:len(self)] = 1
        if os.path.exists(f"{self.outdir}/contigs_assembly.paf.gz"):
            for line in grep(f"{self.name}", f"{self.outdir}/contigs_assembly.paf.gz"):
                alignment = PAF(line)
                if alignment.query_name == self.name:
                    alignments[alignment.query_start:alignment.query_end] = 1
                if alignment.subject_name == self.name:
                    alignments[alignment.subject_start:alignment.subject_end] = 1
        return alignments


    def get_region_depths(self):
        alignments = self.contig_alignments
        regions = alignments.copy()
        regions.split_overlaps()
        region_depths = IntervalTree()
        for region in regions:
            region_depths[region.begin:region.end] = len(alignments[region.begin:region.end])
        
        return sorted(region_depths)


    def get_unique_bases(self):
        unique_bases = len(self)
        for region in self.region_depths:
            if region.data > 1:
                unique_bases -= region.end - region.begin # No need to -1 because end is beyond upper limit
        return unique_bases


    def get_unique_pc(self):
        return self.unique_bases/len(self) * 100


    def assign_ploidy(self, mean, depth):
        ploidy = 0
        prev_diff = None
        while True:
            mean_diff = abs(depth*ploidy - mean)
            if prev_diff is not None and mean_diff > prev_diff:
                ploidy -= 1 # Previous ploidy was a better fit
                break
            prev_diff = mean_diff
            ploidy += 1
        return ploidy


    def get_ploidys(self, median_depth, components=10):

        ploidys = defaultdict(float)
        if len(self.read_depths) < components: # Can't fit model with fewer windows than components
            return ploidys

        model = mixture.BayesianGaussianMixture(n_components=components, max_iter=1000)
        depths = np.array([d.depth for d in self.read_depths]).reshape(-1,1)
        model.fit(depths)

        for i in range(0,len(model.means_)):
            ploidy = self.assign_ploidy(model.means_[i], median_depth/2)
            ploidys[ploidy] += model.weights_[i]

        return ploidys
