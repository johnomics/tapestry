import os, sys, gzip

import logging as log
import networkx as nx

from collections import defaultdict
from intervaltree import Interval, IntervalTree

from .misc import setup_output, PAF
from .misc import minimap2, pigz

def assembly(name):
    nameparts = name.split('_')
    assembly = '_'.join(nameparts[:-1])
    return assembly


class AssemblyReport():
    def __init__(self, assembly, line):
        self.name = assembly

        f = line.rstrip().split('\t')
        self.file = f[0]
        self.length = int(f[1])
        self.gc = float(f[2])
        self.median_read_depth = int(f[3])
        self.unique_bases = int(f[4])
        self.unique_pc = int(f[5])

        self.num_contigs = None


    def __repr__(self):
        report  = f"{self.name:30s}"
        report += f"\t{self.num_contigs}"
        report += f"\t{self.length}"
        report += f"\t{self.gc}"
        report += f"\t{self.median_read_depth}"
        report += f"\t{self.unique_bases}"
        report += f"\t{self.unique_pc}"
        for category in sorted(self.categories, key=lambda x:len(self.categories[x]), reverse=True):
            report += f"\t{category}:{len(self.categories[category])}"
        return report


class ContigReport():
    def __init__(self, line):
        f = line.rstrip().split('\t')
        self.name = f[0]
        self.length = int(f[1])
        self.gc = float(f[2])
        self.median_read_depth = float(f[3])
        self.mean_read_depth = float(f[4])
        self.aligned_primary_pc = float(f[5])
        self.aligned_all_pc = float(f[6])
        self.mean_contig_depth = float(f[7])
        self.tel_start = int(f[8])
        self.tel_end = int(f[9])
        self.mean_start_overhang = int(f[10]) if f[10] != 'None' else None
        self.mean_end_overhang = int(f[11]) if f[11] != 'None' else None
        self.unique_bases = int(f[12])
        self.unique_pc = int(f[13])
        self.category = f[14]
        self.ploidys = '\t'.join(sorted(f[15:]))
    
    def __repr__(self):
        report  = f"{self.name:30s}"
        report += f"\t{self.length}"
        report += f"\t{self.gc}"
        report += f"\t{self.median_read_depth}"
        report += f"\t{self.mean_read_depth}"
        report += f"\t{self.aligned_primary_pc}"
        report += f"\t{self.aligned_all_pc}"
        report += f"\t{self.mean_contig_depth}"
        report += f"\t{self.tel_start}"
        report += f"\t{self.tel_end}"
        report += f"\t{self.mean_start_overhang}"
        report += f"\t{self.mean_end_overhang}"
        report += f"\t{self.unique_bases}"
        report += f"\t{self.unique_pc}"
        report += f"\t{self.category}"
        report += f"\t{self.ploidys}"

        return report


class Stitcher():
    
    def __init__(self, assemblies, outdir, cores):
        self.assemblies = assemblies[0]
        self.outdir = outdir
        self.cores = cores
        self.alignments = defaultdict(lambda: defaultdict(IntervalTree))
        self.contig_graph = nx.Graph()

        setup_output(self.outdir)

        self.align_all()

        self.load_contig_reports()

        self.load_assembly_reports()
        
        self.calculate_assembly_stats()

        self.cluster_contigs()


    def run_minimap2(self, assembly_1, assembly_2, outfilename):
        missing_assembly=False
        for a in assembly_1, assembly_2:
            if not os.path.exists(f"{a}/assembly.fasta"):
                log.error(f"Can't find {a}/assembly.fasta")
                missing_assembly=True
        if missing_assembly:
            sys.exit()

        try:
            log.info(f"Aligning {assembly_1} to {assembly_2}")
            align = (minimap2[f"-xmap-ont", "-2", f"-t{self.cores}", \
                              f"{assembly_1}/assembly.fasta", f"{assembly_2}/assembly.fasta"] | 
                              pigz[f"-p{self.cores}"] > outfilename)
            align()
        except:
            log.error(f"Failed to align {assembly_1} to {assembly_2}")
            sys.exit()


    def load_alignments(self, outfilename):
        alignments = defaultdict(lambda: defaultdict(IntervalTree))
        try:
            with gzip.open(outfilename, 'rt') as paf:
                for line in paf:
                    aln = PAF(line)
                    for name, length in (aln.query_name, aln.query_length), (aln.subject_name, aln.subject_length):
                        if name not in self.contig_graph:
                            self.contig_graph.add_node(name, length=length, assembly=assembly(name))
                        alignments[aln.query_name][aln.subject_name][aln.query_start:aln.query_end] = None
                        alignments[aln.subject_name][aln.query_name][aln.subject_start:aln.subject_end] = None
        except:
            log.error(f"Failed to load alignments from {outfilename}")

        return alignments


    def get_best_hits(self, alignments):
        best_hits = defaultdict(str)
        for contig1 in alignments:
            hit_length = defaultdict(int)
            for contig2 in alignments[contig1]:
                alignments[contig1][contig2].merge_overlaps()
                for a in alignments[contig1][contig2]:
                    hit_length[contig2] += a.end - a.begin
            best_hits[contig1] = sorted(hit_length, key=lambda c: hit_length[c], reverse=True)[0]
        return best_hits


    def fill_contig_graph(self, best_hits):
        for contig in best_hits:
            best_hit_contig = best_hits[contig]
            if best_hits[best_hit_contig] == contig:
                self.contig_graph.add_edge(contig, best_hit_contig)


    def align_pair(self, assembly_1, assembly_2):
        outfilename = f"{self.outdir}/{assembly_1}.{assembly_2}.paf.gz"
        if os.path.exists(outfilename):
            log.info(f"Will use existing {outfilename}")
        else:
            self.run_minimap2(assembly_1, assembly_2, outfilename)

        alignments = self.load_alignments(outfilename)

        best_hits = self.get_best_hits(alignments)

        self.fill_contig_graph(best_hits)


    def align_all(self):
        for a1, assembly_1 in enumerate(sorted(self.assemblies)):
            for a2, assembly_2 in enumerate(sorted(self.assemblies)):
                if a1 >= a2: # Files are the same or already processed
                    continue
                self.align_pair(assembly_1, assembly_2)


    def load_contig_reports(self):
        self.contigs = {}
        for assembly in self.assemblies:
            contig_report = f"{assembly}/contig_report.txt"
            if not os.path.exists(contig_report):
                log.error(f"Can't find {contig_report}")
                sys.exit()
            with open(contig_report) as report:
                for line in report:
                    contig = ContigReport(line)
                    self.contigs[contig.name] = contig


    def load_assembly_reports(self):
        self.assembly_reports = {}
        for assembly in self.assemblies:
            assembly_report = f"{assembly}/assembly_report.txt"
            if not os.path.exists(assembly_report):
                log.error(f"Can't find {assembly_report}")
                sys.exit()
            with open(assembly_report) as report:
                for line in report:
                    assembly_report = AssemblyReport(assembly, line)
                    self.assembly_reports[assembly] = assembly_report


    def calculate_assembly_stats(self):
        for assembly in self.assembly_reports:
            num_contigs = 0
            categories = defaultdict(list)
            for contig in self.contigs:
                if assembly in contig:
                    num_contigs += 1
                    categories[self.contigs[contig].category].append(contig)
            self.assembly_reports[assembly].num_contigs = num_contigs
            self.assembly_reports[assembly].categories = categories


    def cluster_contigs(self):
        clusters = defaultdict(list)
        for c in nx.connected_components(self.contig_graph):
            contigs = self.contig_graph.subgraph(c).nodes()
            clusters[len(contigs)].append(contigs)
        for num_clusters in sorted(clusters, reverse=True):
            print(num_clusters, len(clusters[num_clusters]))
            for cluster in clusters[num_clusters]:
                for contig in sorted(cluster):
                    print(self.contigs[contig])
                print()


    def report(self):
        log.info(f"Generating report")
        with open(f"{self.outdir}/report.txt", 'wt') as report_file:
            for assembly in self.assemblies:
                print(self.assembly_reports[assembly], file=report_file)
        print(f"{self.contig_graph.number_of_nodes()} nodes and {self.contig_graph.number_of_edges()} edges")
