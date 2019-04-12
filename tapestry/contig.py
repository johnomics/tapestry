import os, re, warnings
from statistics import mean
from math import sin, cos
from collections import namedtuple, defaultdict

import numpy as np

import pysam

from intervaltree import Interval, IntervalTree

from sklearn import mixture
from sklearn.exceptions import ConvergenceWarning

from Bio.SeqUtils import GC

from .alignments import Alignments
from .misc import grep, file_exists

# Define process_contig at top level rather than in class so it works with multiprocessing
def process_contig(contig):
    contig.process()
    return contig

def get_ploidy(contig, median_depth=None):
    contig.ploidys = contig.get_ploidys(median_depth)
    contig.ploidy_pc = contig.get_ploidy_pc()
    return contig


class Contig:

    def __init__(self, rec, telomeres, windowsize, outdir, filenames):
        self.name = rec.id
        self.rec = rec
        self.telomeres = telomeres
        self.windowsize = windowsize
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
        report += "\t" + ','.join([f"{p}:{self.ploidy_pc[p]:.2f}" for p in sorted(self.ploidy_pc)])
        report += "\t" + ','.join(self.left_connectors) if self.left_connectors else '\tNone'
        report += "\t" + ','.join(self.right_connectors) if self.right_connectors else '\tNone'
    
        return report


    def json(self):
        return {
            'cluster': self.cluster,
            'longname' : self.name,
            'name': self.name.split('_')[-1], # Remove assembly name
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


    def contig_alignments_json(self):
        plot_row_ends = []
        alignments = []

        for a in sorted(self.contig_alignments, key=lambda a: a.begin):
            assigned_row = None
            for r, row in enumerate(plot_row_ends):
                if row + 1000 < a.begin:
                    assigned_row = r
                    plot_row_ends[r] = a.end
                    break
            if assigned_row is None:
                assigned_row = len(plot_row_ends)
                plot_row_ends.append(a.end)

            alignments.append((a.begin, a.end, assigned_row, a.data))

        return alignments

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
        self.read_alignments = self.plot_read_alignments()

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
        
        max_ploidy = sorted(self.ploidy_pc, key=lambda x:self.ploidy_pc[x], reverse=True)[0] if self.ploidy_pc else 0
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


    def get_read_overhangs(self):

        aligned_length = min(20000, len(self)*0.9)
        start_overhangs = self.alignments.get_start_overhangs(self.name, 1, min(2000, len(self)), aligned_length)
        end_overhangs   = self.alignments.get_end_overhangs(self.name, max(len(self)-2000, 1), len(self), aligned_length)

        mean_start_overhang = int(mean(start_overhangs)) if start_overhangs else None
        mean_end_overhang   = int(mean(end_overhangs)) if end_overhangs else None

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
        for start, end, contig in self.alignments.contig_alignments(self.name):
            alignments[start:end+1] = contig
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


    def get_ploidys(self, median_depth, components=5):

        empty_ploidys = [0] * len(self.read_depths)
        # Can't fit model with fewer windows than components
        if len(self.read_depths) < components: 
            return empty_ploidys

        model = mixture.BayesianGaussianMixture(n_components=components, max_iter=1000)

        depths = np.array(self.read_depths['depth']).reshape(-1,1)

        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        labels = model.fit_predict(depths)
        warnings.resetwarnings()

        haploid_depth = median_depth / 2
        if haploid_depth == 0:
            return empty_ploidys

        ploidys = [int(round((float(model.means_[l]) / haploid_depth))) for l in labels]

        return ploidys


    def get_ploidy_pc(self):
        ploidy_pc = defaultdict(float)
        for p in self.ploidys:
            ploidy_pc[p] += 1/len(self.ploidys)

        return ploidy_pc


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


    def plot_read_alignments(self):
        read_alignments = []
        plot_row_ends = []
        if self.alignments.read_alignments(self.name) is None:
            return read_alignments

        for i, a in self.alignments.read_alignments(self.name).iterrows():
            start_position = a.ref_start - a.left_clip
            end_position = a.ref_end + a.right_clip

            assigned_row = None
            for r, row in enumerate(plot_row_ends):
                if row + 1000 < start_position:
                    assigned_row = r
                    plot_row_ends[r] = end_position
                    break
            if assigned_row is None:
                assigned_row = len(plot_row_ends)
                plot_row_ends.append(end_position)

            # int conversion required because Pandas uses numpy int64, which json doesn't understand
            read_alignments.append([int(x) for x in
                                [start_position,    # read start including left clip
                                 a.ref_start,       # contig alignment start
                                 a.ref_end,         # contig alignment end
                                 end_position,      # read end including right clip
                                 a.mq,              # mapping quality
                                 assigned_row       # y position on plot
                           ]])

        return read_alignments