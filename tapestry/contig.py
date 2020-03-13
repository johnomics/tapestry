# contig.py
# Generate statistics for one contig

# Part of Tapestry
# https://github.com/johnomics/tapestry

# MIT License
# 
# Copyright (c) 2019 John Davey
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import os, re, warnings
from statistics import mean
from collections import Counter, defaultdict

from intervaltree import Interval, IntervalTree

from .alignments import Alignments

# Define process_contig at top level rather than in class so it works with multiprocessing
def process_contig(contig):
    contig.process()
    return contig

def get_ploidy(contig, ploidy_depths=None):
    contig.ploidys = contig.get_ploidys(ploidy_depths)
    contig.ploidy_pc = contig.get_ploidy_pc()
    return contig


class Contig:

    def __init__(self, cid, rec, telomeres, filenames):
        self.id = cid
        self.name = rec.id
        self.rec = rec
        self.telomeres = telomeres
        self.filenames = filenames


    def report(self, assembly_gc):
        report = f"{self.name}"
        report += f"\t{len(self)}"
        report += f"\t{self.gc:.1f}"
        report += f"\t{self.median_read_depth:.1f}"
        report += f"\t{self.tel_start}"
        report += f"\t{self.tel_end}"
        report += f"\t{self.mean_start_overhang}"
        report += f"\t{self.mean_end_overhang}"
        report += f"\t{self.unique_bases}"
        report += f"\t{self.unique_pc:.0f}"
        report += "\t" + ','.join([f"{p}:{self.ploidy_pc[p]:.2f}" for p in sorted(self.ploidy_pc)])

        return report


    def json(self):
        return {
            'id': self.id,
            'group': 'None',
            'name' : self.name,
            'length': len(self),
            'gc': f"{self.gc:.2f}",
            'median_read_depth': int(self.median_read_depth),
            'tel_start': self.tel_start,
            'tel_end': self.tel_end
        }


    def __len__(self):
        return len(self.rec.seq)


    def __lt__(self, other):
        return len(self) < len(other)


    def contig_alignments_json(self):
        plot_row_ends = []
        alignments = []

        for a in sorted(self.contig_alignments, key=lambda a: a.begin):
            contig, contig_start, contig_end = a.data
            if contig == self.name and (a.begin <= contig_end and a.end >= contig_start):
                continue
            alignments.append((a.begin, a.end, self.contig_ids[contig], contig_start, contig_end))

        return alignments

    def process(self):
        # Alignments added here for multithreading
        self.alignments = Alignments(self.filenames['alignments'], self.windowsize)

        self.gc = self.get_gc()
        self.read_depths = self.alignments.depths('read', self.name)
        self.median_read_depth = self.median_depth(self.read_depths)
        self.contig_alignments, self.contig_coverage = self.get_contig_alignments()
        self.mean_start_overhang, self.mean_end_overhang = self.get_read_overhangs()
        self.region_depths = self.get_region_depths()
        self.unique_bases = self.get_unique_bases()
        self.unique_pc = self.get_unique_pc()
        self.tel_start, self.tel_end = self.num_telomeres()
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


    def get_gc(self):
        # Calculate directly because Biopython GC does not account for gaps
        basefreqs = Counter(self.rec.seq)
        for base in "ACGT":
            if base not in basefreqs:
                basefreqs[base] = 0
        gc_bases = sum([basefreqs[b] for b in "GC"])
        acgt_bases = sum([basefreqs[b] for b in "ACGT"])
        gc = 0
        if acgt_bases > 0:
            gc = (sum([basefreqs[b] for b in "GC"]) / sum([basefreqs[b] for b in "ACGT"])) * 100
        return gc


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


    def get_contig_alignments(self):
        alignments = IntervalTree()
        alignments_by_contig = defaultdict(IntervalTree)
        alignments[1:len(self)] = (self.name, 1, len(self))
        for self_start, self_end, contig, contig_start, contig_end in self.alignments.contig_alignments(self.name):
            alignments[self_start:self_end+1] = (contig, contig_start, contig_end)
            alignments_by_contig[contig][self_start:self_end+1] = 1

        coverage = defaultdict(int)
        for contig in alignments_by_contig:
            alignments_by_contig[contig].merge_overlaps()
            coverage[contig] = sum([i.end-i.begin for i in alignments_by_contig[contig]])

        return alignments, coverage


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


    def get_ploidys(self, ploidy_depths):

        empty_ploidys = [0] * len(self.read_depths)
        if sum(ploidy_depths) == 0: 
            return empty_ploidys

        ploidys = []
        for d in self.read_depths['depth']:
            window_ploidy = 0
            min_ploidy_diff = None
            for p, pd in enumerate(ploidy_depths):
                ploidy_diff = abs(d-pd)
                if min_ploidy_diff is None or ploidy_diff < min_ploidy_diff:
                    window_ploidy = p
                    min_ploidy_diff = ploidy_diff
            ploidys.append(window_ploidy)

        return ploidys


    def get_ploidy_pc(self):
        ploidy_pc = defaultdict(float)
        for p in self.ploidys:
            ploidy_pc[p] += 1/len(self.ploidys)

        return ploidy_pc


    def get_neighbour_details(self, neighbour_contig):
        neighbour_id = -1
        neighbour_type = 'L' # Loose
        if neighbour_contig is not None:
            neighbour_id = self.contig_ids[neighbour_contig]
            if neighbour_contig == self.name:
                neighbour_type = 'S' # Self
            else:
                neighbour_type = 'C' # Connection
        return neighbour_id, neighbour_type


    def plot_read_alignments(self):
        read_alignments = []
        
        if not self.readoutput or self.alignments.read_alignments(self.name) is None:
            return read_alignments

        plot_row_ends = []

        for i, a in self.alignments.read_alignments(self.name).iterrows():

            # Only use neighbour distances if contig is the same;
            # if contigs are different, want to show full clips even if the clip aligns elsewhere
            start_distance, end_distance = a.left_clip, a.right_clip
            if a.pre_contig == self.name:
                start_distance = max(a.pre_distance, 0) # Ignore overlapping alignments with negative distances
            if a.post_contig == self.name:
                end_distance = max(a.post_distance, 0)
            
            start_position = a.ref_start - start_distance
            end_position = a.ref_end + end_distance

            assigned_row = None
            for r, row in enumerate(plot_row_ends):
                if row < a.ref_start:
                    assigned_row = r
                    plot_row_ends[r] = a.ref_end
                    break
            if assigned_row is None:
                assigned_row = len(plot_row_ends)
                plot_row_ends.append(end_position)

            pre_contig_id, pre_type = self.get_neighbour_details(a.pre_contig)
            post_contig_id, post_type = self.get_neighbour_details(a.post_contig)

            # int conversion required because Pandas uses numpy int64, which json doesn't understand
            read_alignments.append([int(x) for x in
                                [start_position,    # read start including left clip or pre distance
                                 a.ref_start,       # contig alignment start
                                 a.ref_end,         # contig alignment end
                                 end_position,      # read end including right clip or post distance
                                 a.mq,              # mapping quality
                                 assigned_row+1,    # y position on plot
                                 pre_contig_id,
                                 post_contig_id
                           ]])
            read_alignments[-1].append(pre_type)
            read_alignments[-1].append(post_type)

        return read_alignments
