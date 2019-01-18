import os, sys

import logging as log

from .misc import setup_output

from .misc import minimap2, pigz

class Stitcher():
    
    def __init__(self, assemblies, outdir, cores):
        self.assemblies = assemblies[0]
        self.outdir = outdir
        self.cores = cores

        setup_output(self.outdir)
        
        self.align_all()

    def align_pair(self, assembly_1, assembly_2):
        outfilename = f"{self.outdir}/{assembly_1}.{assembly_2}.paf.gz"
        if os.path.exists(outfilename):
            log.info(f"Will use existing {outfilename}")
        else:
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


    def align_all(self):
        for a1, assembly_1 in enumerate(sorted(self.assemblies)):
            for a2, assembly_2 in enumerate(sorted(self.assemblies)):
                if a1 >= a2: # Files are the same or already processed
                    continue
                self.align_pair(assembly_1, assembly_2)


    def report(self):
        log.info(f"Generating report")
        with open(f"{self.outdir}/report.txt", 'wt') as report_file:
            for assembly in self.assemblies:
                print(f"{assembly}", file=report_file)
