import os, sys, json
import logging as log

from shutil import copyfile
from multiprocessing import Pool
from functools import partial
from statistics import mean, median
from tqdm import tqdm

import networkx as nx

from Bio import SeqIO, motifs
from Bio.Seq import Seq

from jinja2 import Environment, FileSystemLoader
from .assembly_plot import AssemblyPlot
from .contig import Contig, process_contig, get_ploidy

from .misc import flatten, cached_property, setup_output, include_file, report_folder
from .misc import minimap2, samtools, paftools, mosdepth, pigz



class Assembly(AssemblyPlot):

    def __init__(self, assemblyfile, readfile, telomeres, outdir, cores):
        self.assemblyfile = assemblyfile
        self.readfile = readfile
        self.telomeres = [motifs.create([Seq(t[0])]) for t in telomeres] if telomeres else None
        self.outdir = outdir
        self.cores = cores

        setup_output(self.outdir)
        self.contigs = self.load_assembly()

        if self.readfile:
            self.align_to_assembly('reads')
        else:
            log.warning("No read file provided (-r), will skip read metrics unless previous analysis files exist")

        self.align_to_assembly('contigs')

        self.process_contigs()

        self.get_ploidys()

        self.cluster_contigs()

        self.contiglist = sorted(self.contigs, key=lambda c:len(self.contigs[c]), reverse=True)


    def __len__(self):
        return sum([len(self.contigs[c]) for c in self.contigs])


    @cached_property
    def gc(self):
        return mean([float(self.contigs[c].gc) for c in self.contigs])


    @cached_property
    def read_depths(self):
        return flatten([[d.depth for d in self.contigs[c].depths('reads')] for c in self.contigs])


    @cached_property
    def median_depth(self):
        return median(self.read_depths) if self.read_depths else 0


    @cached_property
    def unique_bases(self):
        return sum([self.contigs[c].unique_bases for c in self.contigs])


    @cached_property
    def unique_pc(self):
        return f"{self.unique_bases/len(self) * 100:.0f}"


    def index_assembly(self):
        if os.path.exists(f"{self.outdir}/assembly.fasta.fai"):
           log.info(f"Will use existing {self.outdir}/assembly.fasta.fai")
        else:
            try:
                log.info(f"Indexing assembly")
                samtools("faidx", f"{self.outdir}/assembly.fasta")
            except:
                log.error(f"Can't index assembly!")
                sys.exit()


    def load_assembly(self):
        contigs = {}
        assembly_filename = f"{self.outdir}/assembly.fasta"
        write_assembly = False if os.path.exists(assembly_filename) else True
        if not write_assembly:
            log.info(f"Will use existing {assembly_filename}")

        try:
            log.info(f"Loading genome assembly")
            assembly_out = open(assembly_filename, 'w') if write_assembly else None
            for rec in SeqIO.parse(open(self.assemblyfile, 'r'), "fasta"):
                rec.seq = rec.seq.upper()
                rec.id = f"{self.outdir}_{rec.id}"
                contigs[rec.id] = Contig(rec, self.telomeres, self.outdir)
                if write_assembly:
                    SeqIO.write(rec, assembly_out, "fasta")
        except IOError:
            log.error(f"Can't load assembly from file {self.assemblyfile}!")
            sys.exit()

        self.index_assembly()

        return contigs


    def make_bam(self, aligntype):
        if os.path.exists(f"{self.outdir}/{aligntype}_assembly.bam"):
            log.info(f"Will use existing {self.outdir}/{aligntype}_assembly.bam")
        else:
            inputfile = x_option = None
            if aligntype == 'reads':
                inputfile = os.path.abspath(self.readfile)
                x_option = 'map-ont'
            elif aligntype == 'contigs':
                inputfile = f"{self.outdir}/assembly.fasta"
                x_option = 'ava-ont'
            else:
                log.error(f"Don't know how to align {aligntype}")
                sys.exit()
            log.info(f"Aligning {aligntype} {inputfile} to assembly")

            # samtools uses cores-1 because -@ specifies additional cores and defaults to 0
            try:
                align = minimap2[f"-x{x_option}", "-a", "-2", f"-t{self.cores}", \
                                       f"{self.outdir}/assembly.fasta", inputfile] | \
                              samtools["sort", f"-@{self.cores-1}", \
                                       f"-o{self.outdir}/{aligntype}_assembly.bam"]
                align()
            except:
                log.error(f"Failed to align {inputfile} to {self.outdir}/assembly.fasta")
                sys.exit()


    def index_bam(self, aligntype):
        if os.path.exists(f"{self.outdir}/{aligntype}_assembly.bam.bai"):
            log.info(f"Will use existing {self.outdir}/{aligntype}_assembly.bam.bai")
        else:
            try:
                samtools("index", f"-@{self.cores-1}", f"{self.outdir}/{aligntype}_assembly.bam")
            except:
                log.error(f"Failed to index {self.outdir}/{aligntype}_assembly.bam")


    def make_paf(self, aligntype):
        if os.path.exists(f"{self.outdir}/{aligntype}_assembly.paf.gz"):
            log.info(f"Will use existing {self.outdir}/{aligntype}_assembly.paf.gz")
        else:
            try:
                samtopaf = (samtools["view", "-h", f"{self.outdir}/{aligntype}_assembly.bam"] | \
                           paftools["sam2paf", "-"] | pigz[f"-p{self.cores}"] > f"{self.outdir}/{aligntype}_assembly.paf.gz")
                samtopaf()
            except:
                log.error(f"Failed to convert {self.outdir}/{aligntype}_assembly.bam to {self.outdir}/{aligntype}_assembly.paf.gz")
                sys.exit()


    def run_mosdepth(self, filestub):
        if os.path.exists(f"{self.outdir}/{filestub}.regions.bed.gz"):
            log.info(f"Will use existing {filestub} mosdepth output")
            return

        log.info(f"Running mosdepth for {filestub}")
        
        try:
            mosdepth("-b1000", "-n", \
                     f"{self.outdir}/{filestub}", \
                     f"{self.outdir}/{filestub}.bam")
        except:
            log.error(f"Failed to run mosdepth for {self.outdir}/{filestub}.bam")
            sys.exit()


    def align_to_assembly(self, aligntype):
        self.make_bam(aligntype)
        self.index_bam(aligntype)
        self.run_mosdepth(f"{aligntype}_assembly")
        if aligntype=="contigs":
            self.make_paf(aligntype)


    def process_contigs(self):
        log.info(f"Processing {len(self.contigs)} contigs")
        with Pool(self.cores) as p:
            for contig in tqdm(p.imap(process_contig, self.contigs.values()), total=len(self.contigs)):
                self.contigs[contig.name] = contig


    def get_ploidys(self):
        log.info(f"Calculating ploidy estimates")
        fit_ploidy = partial(get_ploidy, median_depth=self.median_depth)
        with Pool(self.cores) as p:
            for contig in tqdm(p.imap(fit_ploidy, self.contigs.values()), total=len(self.contigs)):
                self.contigs[contig.name] = contig


    def cluster_contigs(self):
        log.info(f"Clustering contigs")
        self.assembly_graph = nx.Graph()
        for contig in self.contigs:
            self.assembly_graph.add_node(contig)

        for contig in self.contigs:
            for connector in self.contigs[contig].left_connectors + self.contigs[contig].right_connectors:
                connector_contig = connector.split(':')[0]

                if self.assembly_graph.has_edge(contig, connector_contig):
                    self.assembly_graph[contig][connector_contig]['links'] += 1
                else:
                    self.assembly_graph.add_edge(contig, connector_contig, links=1)

        self.components = []
        for c in nx.connected_components(self.assembly_graph):
            self.components.append(self.assembly_graph.subgraph(c).nodes())

        for i, component in enumerate(sorted(self.components, key=lambda c: len(c), reverse=True)):
            for contig in sorted(component):
                self.contigs[contig].cluster = i+1


    def contig_report(self):
        log.info(f"Generating contig report")

        try:
            with open(f"{self.outdir}/contig_report.txt", 'wt') as report_file, \
                 open(f"{self.outdir}/redundancy.txt", 'wt') as redundancy_file:
                for contigname in sorted(self.contigs, key=lambda c: (self.contigs[c].cluster, -len(self.contigs[c]))):
                    print(self.contigs[contigname].report(self.gc), file=report_file)
                    print(self.contigs[contigname].redundancy_report(), file=redundancy_file)
        except IOError:
            log.error(f"Could not write contig reports")


    def assembly_report(self):
        log.info(f"Generating assembly report")
        with open(f"{self.outdir}/assembly_report.txt", 'wt') as report_file:
            print(f"{self.assemblyfile}\t{len(self)}\t{self.gc:.1f}\t{self.median_depth:.0f}\t{self.unique_bases}\t{self.unique_pc}", file=report_file)


    def html_report(self):
        env = Environment(loader=FileSystemLoader(report_folder))
        env.globals['include_file'] = include_file
        template = env.get_template('template.html')
        with open(f"{self.outdir}/tapestry_report.html", 'wt') as html_report:
            print(template.render(contigs=json.dumps([self.contigs[c].json() for c in self.contigs])), file=html_report)


    def report(self):
        log.info(f"Generating reports")
        self.contig_report()
        self.assembly_report()
        self.html_report()


    def plot(self):
        log.info(f"Plotting")
        self.lengthplot()
        if self.readfile:
            self.depthplot_full()
            self.depthplot_summary()
