import os, sys, errno
import logging as log
from shutil import copyfile
from multiprocessing import Pool
from statistics import mean

from Bio import SeqIO, motifs
from Bio.Seq import Seq

from .assembly_plot import AssemblyPlot
from .contig import Contig, process_contig

from .misc import flatten, cached_property
from .misc import minimap2, samtools, paftools, mosdepth, pigz



class Assembly(AssemblyPlot):

    def __init__(self, assemblyfile, readfile, telomeres, outdir, cores):
        self.assemblyfile = assemblyfile
        self.readfile = readfile
        self.telomeres = [motifs.create([Seq(t[0])]) for t in telomeres] if telomeres else None
        self.outdir = outdir
        self.cores = cores

        self.setup_output()
        self.prepare_genome()
        self.contigs = self.load_genome()

        if self.readfile:
            self.align_to_assembly('reads')
        else:
            log.warning("No read file provided (-r), will skip read metrics unless previous analysis files exist")

        self.align_to_assembly('contigs')

        self.process_contigs()
        
        self.contiglist = sorted(self.contigs, key=lambda c:len(self.contigs[c]), reverse=True)


    def __len__(self):
        return sum([len(self.contigs[c]) for c in self.contigs])


    @cached_property
    def gc(self):
        return mean([float(self.contigs[c].gc) for c in self.contigs])


    @cached_property
    def median_depth(self):
        depths = flatten([[d.depth for d in self.contigs[c].depths('reads')] for c in self.contigs])
        return depths[int(len(depths)/2)] if depths else 0

    @cached_property
    def unique_bases(self):
        return sum([self.contigs[c].unique_bases for c in self.contigs])


    @cached_property
    def unique_pc(self):
        return f"{self.unique_bases/len(self) * 100:.0f}"


    def setup_output(self):
        try:
            os.mkdir(self.outdir)
            log.info(f"Created output directory {self.outdir}")
        except OSError as exc:
            if exc.errno == errno.EEXIST:
                log.warning(f"Output directory {self.outdir} found, will use existing analysis files if present, but overwrite reports")
            else:
                raise


    def prepare_genome(self):
        if os.path.exists(f"{self.outdir}/assembly.fasta"):
            log.info(f"Will use existing {self.outdir}/assembly.fasta")
        else:
            try:
                log.info(f"Copying assembly to {self.outdir}")
                copyfile(self.assemblyfile, f"{self.outdir}/assembly.fasta") 
            except:
                 log.error(f"Can't copy assembly to {self.outdir}")
                 sys.exit()

        if os.path.exists(f"{self.outdir}/assembly.fasta.fai"):
           log.info(f"Will use existing {self.outdir}/assembly.fasta.fai")
        else:
            try:
                log.info(f"Indexing assembly")
                samtools("faidx", f"{self.outdir}/assembly.fasta")
            except:
                log.error(f"Can't index assembly!")
                sys.exit()


    def load_genome(self):
        contigs = {}
        try:
            log.info(f"Loading genome assembly")
            for rec in SeqIO.parse(open(f"{self.outdir}/assembly.fasta", 'r'), "fasta"):
                rec.seq = rec.seq.upper()
                contigs[rec.id] = Contig(rec, self.telomeres, self.outdir)
        except IOError:
            log.error(f"Can't load assembly from file {self.assemblyfile}!")
            sys.exit()

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
                samtools("index", f"-@{self.cores-1}", f"{self.outdir}/{aligntype}_assembly.bam")
            except:
                log.error(f"Failed to align {inputfile} to {self.outdir}/assembly.fasta")
                sys.exit()


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
            mosdepth(f"-t{self.cores}", "-b1000", "-n", \
                     f"{self.outdir}/{filestub}", \
                     f"{self.outdir}/{filestub}.bam")
        except:
            log.error(f"Failed to run mosdepth for {self.outdir}/{filestub}.bam")
            sys.exit()


    def align_to_assembly(self, aligntype):
        self.make_bam(aligntype)
        self.run_mosdepth(f"{aligntype}_assembly")
        if aligntype=="contigs":
            self.make_paf(aligntype)


    def process_contigs(self):
        log.info(f"Processing contigs")
        with Pool(self.cores) as p:
            for contig in p.map(process_contig, self.contigs.values()):
                self.contigs[contig.name] = contig


    def contig_report(self):
        log.info(f"Generating contig report")

        try:
            with open(f"{self.outdir}/contig_report.txt", 'wt') as report_file, \
                 open(f"{self.outdir}/redundancy.txt", 'wt') as redundancy_file:
                for contigname in self.contigs:
                    print(self.contigs[contigname], file=report_file)
                    print(self.contigs[contigname].redundancy_report(), file=redundancy_file)
        except IOError:
            log.error(f"Could not write contig reports")


    def assembly_report(self):
        log.info(f"Generating assembly report")
        with open(f"{self.outdir}/assembly_report.txt", 'wt') as report_file:
            print(f"{self.assemblyfile}\t{len(self)}\t{self.gc:.1f}\t{self.median_depth:.0f}\t{self.unique_bases}\t{self.unique_pc}", file=report_file)


    def report(self):
        log.info(f"Generating reports")
        self.contig_report()
        self.assembly_report()


    def plot(self):
        log.info(f"Plotting")
        self.lengthplot()
        if self.readfile:
            self.depthplot_full()
            self.depthplot_summary()
