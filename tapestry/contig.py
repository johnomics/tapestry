import os, re, warnings
from collections import namedtuple, defaultdict

import numpy as np

import pysam

from intervaltree import Interval, IntervalTree

from sklearn import mixture
from sklearn.exceptions import ConvergenceWarning

from Bio.SeqUtils import GC

from .alignments import Alignments
from .misc import grep, PAF, file_exists

# Define process_contig at top level rather than in class so it works with multiprocessing
def process_contig(contig):
    contig.process()
    return contig

def get_ploidy(contig, median_depth=None):
    contig.ploidys = contig.get_ploidys(median_depth)
    return contig


class Contig:

    def __init__(self, rec, telomeres, outdir, filenames):
        self.name = rec.id
        self.rec = rec
        self.telomeres = telomeres
        self.outdir = outdir
        self.filenames = filenames


    def report(self, assembly_gc):
        report = f"{self.cluster}"
        report += f"\t{self.name}"
        report += f"\t{len(self)}"
        report += f"\t{self.gc:.1f}"
        report += f"\t{self.median_read_depth:.1f}"
        report += f"\t{self.mean_read_depth:.1f}"
        report += f"\t{self.aligned_primary_pc:.1f}"
        report += f"\t{self.aligned_all_pc:.1f}"
        report += f"\t{self.mean_contig_depth:.1f}"
        report += f"\t{self.tel_start}"
        report += f"\t{self.tel_end}"
        report += f"\t{self.mean_start_overhang}"
        report += f"\t{self.mean_end_overhang}"
        report += f"\t{self.unique_bases}"
        report += f"\t{self.unique_pc:.0f}"
        report += f"\t{self.category(assembly_gc)}"
        report += "\t" + ','.join([f"{p}:{self.ploidys[p]:.2f}" for p in sorted(self.ploidys)])
        report += "\t" + ','.join(self.left_connectors) if self.left_connectors else '\tNone'
        report += "\t" + ','.join(self.right_connectors) if self.right_connectors else '\tNone'
    
        return report


    def json(self):
        return {
            'cluster': self.cluster,
            'name': self.name,
            'length': len(self),
            'gc': self.gc,
            'tel_start': self.tel_start,
            'tel_end': self.tel_end
        }


    def __len__(self):
        return len(self.rec.seq)


    def __lt__(self, other):
        return len(self) < len(other)


    def redundancy_report(self):
        regions = ""
        for region in self.region_depths:
            regions += f"{self.name}\t{region.begin}\t{region.end}\t{region.end-region.begin}\t{region.data}\n"
        return regions.rstrip()


    def process(self):
        # Alignments added here for multithreading
        self.alignments = Alignments(self.filenames['alignments'])

        self.gc = self.get_gc()
        self.read_depths = self.alignments.depths('read', self.name)
        self.contig_depths = self.alignments.depths('contig', self.name)
        self.mean_contig_depth = self.mean_depth(self.contig_depths)
        self.mean_read_depth = self.mean_depth(self.read_depths)
        self.median_read_depth = self.median_depth(self.read_depths)
        self.contig_alignments = self.get_contig_alignments()
        self.mean_start_overhang, self.mean_end_overhang = self.get_read_overhangs()
        self.region_depths = self.get_region_depths()
        self.unique_bases = self.get_unique_bases()
        self.unique_pc = self.get_unique_pc()
        self.tel_start, self.tel_end = self.num_telomeres()
        self.left_connectors, self.right_connectors = self.get_connectors()
        self.aligned_primary_pc, self.aligned_all_pc = self.get_read_alignment_stats()
#        self.assess_strength(5000, 2500)

        # Alignments work is done; they cannot be pickled, so clean up before return
        del(self.alignments)

    def completeness(self):
        completeness = ''
        if self.tel_start > 0 and self.mean_start_overhang is not None and self.mean_start_overhang < 250:
            completeness += 'L'
        if self.tel_end   > 0 and self.mean_end_overhang is not None and self.mean_end_overhang   < 250:
            completeness += 'R'
        if completeness == 'LR':
            completeness = 'C'
        return completeness if completeness else '-'


    def category(self, assembly_gc):
        category=''
        if abs(self.gc - assembly_gc) < 2:
            category += 'N'
        else:
            category += '-'

        category += self.completeness()
        
        max_ploidy = sorted(self.ploidys, key=lambda x:self.ploidys[x], reverse=True)[0] if self.ploidys else 0
        if max_ploidy > 4:
            max_ploidy = 'R'
        category += f"{max_ploidy}"

        return category


    def get_gc(self):
        return GC(self.rec.seq)


    def mean_depth(self, depths):
        return depths['depth'].mean() if depths is not None else 0


    def median_depth(self, depths):
        return depths['depth'].median() if depths is not None else 0


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
                    bam, 0, min(1000, len(self)),
                    lambda aln: aln.query_alignment_start-aln.reference_start)
            mean_end_overhang = self.get_read_overhang(
                    bam, max(len(self)-1001,0), len(self)-1,
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


    def get_read_alignment_stats(self):
        
        aligned_primary_bases = all_primary_bases = aligned_non_primary_bases = 0

        count_bases = self.alignments.get_contig_read_counts(self.name)
        if count_bases is not None:
            aligned_non_primary_bases = count_bases.loc['secondary','aligned_length'] + count_bases.loc['supplementary', 'aligned_length']
            all_primary_bases = count_bases.loc['primary', 'read_length']
            aligned_primary_bases = count_bases.loc['primary', 'aligned_length']

        all_aligned_bases = aligned_primary_bases + aligned_non_primary_bases
        
        aligned_primary_pc = aligned_all_pc = 0
        if all_primary_bases > 0:
            aligned_primary_pc = aligned_primary_bases/all_primary_bases*100
        if all_aligned_bases > 0:
            aligned_all_pc = aligned_primary_bases/all_aligned_bases*100

        return aligned_primary_pc, aligned_all_pc


    def get_contig_alignments(self):
        alignments = IntervalTree()
        alignments[1:len(self)] = 1
        for start, end in self.alignments.get_contig_alignments(self.name):
            alignments[start:end+1] = 1
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
        depths = np.array(self.read_depths['depth']).reshape(-1,1)
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        model.fit(depths)
        warnings.resetwarnings()

        for i in range(0,len(model.means_)):
            ploidy = self.assign_ploidy(model.means_[i], median_depth/2)
            ploidys[ploidy] += model.weights_[i]

        return ploidys


    def get_aln_length(self, cigar):
        # Just gets the Matching bases, which is not accurate but good enough for finding connectors
        aln_length = 0
        for cigar_op in re.findall('\d+\w', cigar):
            if cigar_op.endswith('M'):
                aln_length += int(cigar_op[:-1])
        return aln_length


    def get_connections(self, bam, region_start, region_end):
        connections = defaultdict(IntervalTree)
        region_reads = defaultdict(int)

        if region_start < 0:
            region_start = 0

        if region_end > len(self):
            region_end = len(self)

        for aln in bam.fetch(self.name, region_start, region_end):
            region_reads[aln.query_name] = 1
            if aln.has_tag('SA'):
                aln_dir = '-' if aln.is_reverse else '+'
                aln_connections = aln.get_tag('SA').split(';')[:-1]
                for c in aln_connections:
                    contig, start, direction, cigar, *c_args = c.split(',')
                    if contig == self.name:
                        continue
                    aln_length = self.get_aln_length(cigar)
                    connections[contig][int(start):int(start)+aln_length] = 1

        return connections, region_reads


    def get_region_connectors(self, bam, region_start, region_end):
        connectors = []

        connections, region_reads = self.get_connections(bam, region_start, region_end)

        for contig in connections:
            connections[contig].merge_overlaps()

            contig_reads = defaultdict(int)
            for interval in connections[contig]:
                for interval_aln in bam.fetch(contig, interval.begin, interval.end):
                    contig_reads[interval_aln.query_name] = 1

            connecting_reads = set(region_reads).intersection(set(contig_reads))

            region_connecting_reads_pc = len(connecting_reads)/len(region_reads) if region_reads else 0
            contig_connecting_reads_pc = len(connecting_reads)/len(contig_reads) if contig_reads else 0
            if region_connecting_reads_pc >= 0.65 and contig_connecting_reads_pc >= 0.1:
                connectors.append(f"{contig}:{region_connecting_reads_pc:.2f}:{contig_connecting_reads_pc:.2f}")

        return connectors


    def get_connectors(self):
        left_connectors = right_connectors = []

        if file_exists(self.filenames['reads_bam']):
            bam = pysam.AlignmentFile(self.filenames['reads_bam'], 'rb')

            if self.completeness() in ['R', '-']:
                left_connectors = self.get_region_connectors(bam, 0, 10000)
            if self.completeness() in ['L', '-']:
                right_connectors = self.get_region_connectors(bam, len(self)-10000, len(self)-1)

        return left_connectors, right_connectors

#    def assess_strength(self, region_width, region_step):
#        conn, db = load_reads_database(self.filenames['reads_db'])
#        regions = []
#        for region_start in range(0, len(self), region_step):
#            region_end = min(len(self), region_start + region_width)
#            through_reads, through_av_len, left_reads, left_av_clip, right_reads, right_av_clip = \
#                get_reads_from_region(conn, db, self.name, region_start, region_end)
#            regions.append([self.name, region_start, region_end, through_reads, through_av_len, left_reads, left_av_clip, right_reads, right_av_clip])
#            if region_end == len(self):
#                break
#        conn.close()
