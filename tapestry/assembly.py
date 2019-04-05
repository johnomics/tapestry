import os, sys, json
import logging as log
import networkx as nx

from multiprocessing import Pool
from functools import partial
from statistics import mean, median
from gzip import open as gzopen

from Bio import SeqIO, motifs
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from jinja2 import Environment, FileSystemLoader
from .assembly_plot import AssemblyPlot

from .contig import Contig, process_contig, get_ploidy
from .db import build_reads_database
from .misc import flatten, cached_property, setup_output, include_file, report_folder, tapestry_tqdm, file_exists
from .misc import minimap2, samtools, paftools, mosdepth, pigz

filenames = {
    'assembly'         : 'assembly.fasta', 
    'assembly_index'   : 'assembly.fasta.fai',
    'sampled_reads'    : 'reads.fastq.gz', 
    'reads_bam'        : 'reads_assembly.bam',
    'reads_index'      : 'reads_assembly.bam.bai',
    'reads_mosdepth'   : 'reads_assembly.regions.bed.gz',
    'reads_db'         : 'aligned_reads.db',
    'contigs_bam'      : 'contigs_assembly.bam',
    'contigs_index'    : 'contigs_assembly.bam.bai',
    'contigs_paf'      : 'contigs_assembly.paf.gz',
    'contigs_mosdepth' : 'contigs_assembly.regions.bed.gz'
}



class Assembly(AssemblyPlot):

    def __init__(self, assemblyfile, readfile, telomeres, outdir, cores, depth, minreadlength):
        self.assemblyfile = assemblyfile
        self.readfile = readfile
        self.telomeres = [motifs.create([Seq(t[0])]) for t in telomeres] if telomeres else None
        self.outdir = outdir
        self.cores = cores
        self.depth = depth
        self.minreadlength = minreadlength

        setup_output(self.outdir)
        self.filenames = {file_key:f"{self.outdir}/{filenames[file_key]}" for file_key in filenames}

        self.contigs = self.load_assembly()

        if self.readfile:
            self.sample_reads()
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


    def load_assembly(self):
        contigs = {}
        no_assembly_found = not file_exists(self.filenames['assembly'])
        if no_assembly_found:
            log.info(f"Will use existing {self.filenames['assembly']}")

        try:
            log.info(f"Loading genome assembly")
            assembly_out = open(self.filenames['assembly'], 'w') if no_assembly_found else None
            for rec in SeqIO.parse(open(self.assemblyfile, 'r'), "fasta"):
                rec.seq = rec.seq.upper()
                rec.id = f"{self.outdir}_{rec.id}"
                contigs[rec.id] = Contig(rec, self.telomeres, self.outdir, self.filenames)
                if no_assembly_found:
                    SeqIO.write(rec, assembly_out, "fasta")
        except IOError:
            log.error(f"Can't load assembly from file {self.assemblyfile}!")
            sys.exit()

        self.index_assembly()

        return contigs


    def index_assembly(self):
        if file_exists(self.filenames['assembly_index'], deps=[self.filenames['assembly']]):
           log.info(f"Will use existing {self.filenames['assembly_index']}")
        else:
            try:
                log.info(f"Indexing assembly")
                samtools("faidx", self.filenames['assembly'])
            except:
                log.error(f"Can't index assembly {self.filenames['assembly']}!")
                sys.exit()


    def sample_reads(self):
        if file_exists(self.filenames['sampled_reads']):
            log.info(f"Will use existing {self.filenames['sampled_reads']}")
        else:
            if self.depth == 0:
                os.symlink(os.path.abspath(self.readfile), self.filenames['sampled_reads'])
                log.info(f"Using all reads as depth option is {self.depth}")
            else:
                log.info(f"Sampling {self.depth} times coverage of {len(self)/1000000:.1f} Mb assembly from >{self.minreadlength}bp reads in {self.readfile}")
                bases = len(self) * self.depth
                with gzopen(self.readfile, 'rt') as all_reads, gzopen(self.filenames['sampled_reads'], "wt") as sampled_reads:
                    sampled_bases = 0
                    read_count = 0
                    readset = []
                    for read in SeqIO.parse(all_reads, "fastq"):
                        if len(read.seq) < self.minreadlength:
                            continue
                        readset.append(read)
                        read_count += 1
                        if read_count % 10000 == 0:
                            SeqIO.write(readset, sampled_reads, "fastq")
                            readset = []
                        sampled_bases += len(read.seq)
                        if sampled_bases > bases:
                            break
                    SeqIO.write(readset, sampled_reads, "fastq")
                log.info(f"Wrote {read_count} reads ({sampled_bases} bases) to {self.filenames['sampled_reads']}")


    def make_bam(self, aligntype):
        bam_filename = self.filenames[f'{aligntype}_bam']
        if file_exists(bam_filename, deps=[self.filenames['sampled_reads'], self.filenames['assembly']]):
            log.info(f"Will use existing {bam_filename}")
        else:
            inputfile = x_option = None
            if aligntype == 'reads':
                inputfile = self.filenames['sampled_reads']
                x_option = 'map-ont'
            elif aligntype == 'contigs':
                inputfile = self.filenames['assembly']
                x_option = 'ava-ont'
            else:
                log.error(f"Don't know how to align {aligntype}")
                sys.exit()
            log.info(f"Aligning {aligntype} {inputfile} to assembly")

            # samtools uses cores-1 because -@ specifies additional cores and defaults to 0
            try:
                align = minimap2[f"-x{x_option}", "-a", "-2", f"-t{self.cores}", \
                                       self.filenames['assembly'], inputfile] | \
                              samtools["sort", f"-@{self.cores-1}", \
                                       f"-o{bam_filename}"]
                align()
            except:
                log.error(f"Failed to align {inputfile} to {self.filenames['assembly']}")
                sys.exit()


    def index_bam(self, aligntype):
        bam_filename = self.filenames[f'{aligntype}_bam']
        bai_filename = self.filenames[f'{aligntype}_index']
        if file_exists(bai_filename, deps=[bam_filename]):
            log.info(f"Will use existing {bai_filename}")
        else:
            try:
                log.info(f"Indexing {bam_filename}")
                samtools("index", f"-@{self.cores-1}", bam_filename)
            except:
                log.error(f"Failed to index {bam_filename}")


    def make_paf(self, aligntype):
        paf_filename = self.filenames[f"{aligntype}_paf"]
        bam_filename = self.filenames[f"{aligntype}_bam"]
        if file_exists(paf_filename, deps=[bam_filename]):
            log.info(f"Will use existing {paf_filename}")
        else:
            try:
                log.info(f"Converting {bam_filename} to {paf_filename}")
                samtopaf = (samtools["view", "-h", bam_filename] | \
                           paftools["sam2paf", "-"] | pigz[f"-p{self.cores}"] > paf_filename)
                samtopaf()
            except:
                log.error(f"Failed to convert {bam_filename} to {paf_filename}")
                sys.exit()


    def run_mosdepth(self, aligntype):
        mos_filename = self.filenames[f'{aligntype}_mosdepth']
        bam_filename = self.filenames[f'{aligntype}_bam']
        if file_exists(mos_filename, deps=[bam_filename]):
            log.info(f"Will use existing {mos_filename} mosdepth output")
            return

        try:
            log.info(f"Running mosdepth for {aligntype}")
            mosdepth("-b1000", "-n", \
                     f"{self.outdir}/{aligntype}_assembly", \
                     bam_filename)
        except:
            log.error(f"Failed to run mosdepth for {bam_filename}")
            sys.exit()


    def align_to_assembly(self, aligntype):
        self.make_bam(aligntype)
        self.index_bam(aligntype)
        self.run_mosdepth(aligntype)
        if aligntype=="contigs":
            self.make_paf(aligntype)
        elif aligntype=="reads":
            build_reads_database(self.filenames['reads_bam'], self.filenames['reads_db'], self.contigs)


    def process_contigs(self):
        log.info(f"Processing {len(self.contigs)} contigs")
        with Pool(self.cores) as p:
            for contig in tapestry_tqdm(p.imap(process_contig, self.contigs.values()), total=len(self.contigs), desc="Processing contigs"):
                self.contigs[contig.name] = contig


    def get_ploidys(self):
        log.info(f"Calculating ploidy estimates")
        fit_ploidy = partial(get_ploidy, median_depth=self.median_depth)
        with Pool(self.cores) as p:
            for contig in tapestry_tqdm(p.imap(fit_ploidy, self.contigs.values()), total=len(self.contigs), desc="Ploidy estimates"):
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
            print(template.render(
                    contigs=json.dumps([self.contigs[c].json() for c in self.contigs]),
                    threads=json.dumps([self.contigs[c].json() for c in self.contigs])
                 ),
                 file=html_report)


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
