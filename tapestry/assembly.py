# assembly.py
# Generate statistics and reports for an assembly

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


import os, sys, json, datetime
import logging as log

from multiprocessing import Pool
from functools import partial
from statistics import mean, median
from collections import defaultdict
from gzip import open as gzopen
from tqdm import tqdm
from math import log10

from Bio import SeqIO, motifs
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from jinja2 import Environment, FileSystemLoader

from .contig import Contig, process_contig, get_ploidy
from .alignments import Alignments
from .misc import cached_property, setup_output, include_file, report_folder, tapestry_tqdm, file_exists, is_gz_file
from .misc import minimap2, samtools
from ._version import __version__

filenames = {
    'assembly'         : 'assembly.fasta', 
    'sampled_reads'    : 'reads.fastq.gz', 
    'reads_bam'        : 'reads_assembly.bam',
    'reads_index'      : 'reads_assembly.bam.bai',
    'contigs_bam'      : 'contigs_assembly.bam',
    'contigs_index'    : 'contigs_assembly.bam.bai',
    'alignments'       : 'alignments.db'
}



class Assembly():

    def __init__(self, assemblyfile, readfile, telomeres, outdir, cores, coverage, minreadlength, windowsize, forcereadoutput, mincontigalignment):
        self.assemblyfile = assemblyfile
        self.readfile = readfile
        self.telomere_seqs = ' '.join(telomeres[0]) if telomeres else ''
        self.telomeres = [motifs.create([Seq(t[0])]) for t in telomeres] if telomeres else None
        self.outdir = outdir
        self.basedir = os.path.basename(outdir)
        self.cores = cores
        self.coverage = coverage
        self.minreadlength = minreadlength
        self.forcereadoutput = forcereadoutput
        self.min_contig_alignment = mincontigalignment

        setup_output(self.outdir)
        self.filenames = {file_key:f"{self.outdir}/{filenames[file_key]}" for file_key in filenames}

        self.contigs = self.load_assembly()
        
        self.windowsize = self.get_windowsize(windowsize)
        self.readoutput = len(self) < 50000000 or self.forcereadoutput

        for contig in self.contigs:
            self.contigs[contig].readoutput = self.readoutput
            self.contigs[contig].windowsize = self.windowsize
            
        if not self.readoutput or len(self.contigs) > 300:
            log.warning("Tapestry is designed for small genome assemblies (<50 Mb) that are close to complete (<300 contigs).")
            log.warning(f"It may not perform well for your {len(self)/1000000:.1f} Mb assembly with {len(self.contigs)} contigs.")
            log.warning("The run may take a long time, and the report may be too large to load in your web browser.")
            log.warning("To save time and space, no read alignments will be output for this assembly (use -f to force read alignment output).")
            log.warning("To improve performance, consider increasing window size (-w) and minimum contig alignment (-m), and decreasing read depth (-d).")

        if self.readfile:
            self.sample_reads()
            self.align_to_assembly('reads')
        else:
            log.warning("No read file provided (-r), will skip read metrics unless previous analysis files exist")

        self.align_to_assembly('contigs')

        self.alignments = self.load_alignments()

        self.process_contigs()

        self.get_ploidys()


    def __len__(self):
        return sum([len(self.contigs[c]) for c in self.contigs])


    @cached_property
    def gc(self):
        return mean([float(self.contigs[c].gc) for c in self.contigs])


    @cached_property
    def n50_length(self):
        n50_length = 0
        n50_cumul = 0
        for c in sorted(self.contigs, key=lambda x: -len(self.contigs[x])):
            n50_length = len(self.contigs[c])
            n50_cumul += n50_length
            if n50_cumul >= len(self)/2:
                break
        return n50_length

    @cached_property
    def read_depths(self):
        return self.alignments.depths('read')['depth']


    @cached_property
    def median_depth(self):
        return self.read_depths.median() if self.read_depths is not None else 0

    @cached_property
    def ploidy_depths(self):
        return [d*self.median_depth for d in [0, 0.5, 1, 1.5, 2, 2.5]]

    def get_windowsize(self, windowsize=None):
        if windowsize is None:
            order = 10 ** int(log10(self.n50_length+1))
            multiplier = int(self.n50_length/order)
            windowsize = int((order * multiplier) / 20)
            if windowsize < 10000:
                windowsize = 10000
            log.info(f"Ploidy window size set to {windowsize} based on assembly N50 length {self.n50_length}")
        return windowsize

    def options(self):
        return [
            {'option': 'Tapestry version',         'value': __version__        },
            {'option': 'Report generation time',   'value': datetime.datetime.now().strftime('%d %B %Y %H:%M:%S') },
            {'option': 'Assembly file',            'value': os.path.basename(self.assemblyfile)  },
            {'option': 'Reads file',               'value': os.path.basename(self.readfile)      },
            {'option': 'Telomeres',                'value': self.telomere_seqs },
            {'option': 'Output directory',         'value': self.outdir        },
            {'option': 'Genome coverage',          'value': self.coverage      },
            {'option': 'Minimum read length',      'value': self.minreadlength },
            {'option': 'Window size',              'value': self.windowsize    },
            {'option': 'Minimum contig alignment', 'value': self.min_contig_alignment },
            {'option': 'Read alignments included', 'value': self.readoutput  }
        ]

    def load_assembly(self):
        contigs = {}
        contig_ids = {}
        
        try:
            log.info("Loading genome assembly")
            if is_gz_file(self.assemblyfile):
                self.filenames['assembly'] = f"{self.outdir}/assembly.fasta.gz"

            if file_exists(self.filenames['assembly']):
                log.info(f"Found {self.filenames['assembly']}, will not overwrite")
            else:
                os.symlink(os.path.abspath(self.assemblyfile), self.filenames['assembly'])

            contig_id = 0
            
            assembly_in = None
            if is_gz_file(self.filenames['assembly']):
                assembly_in = gzopen(self.filenames['assembly'], 'rt')
            else:
                assembly_in = open(self.filenames['assembly'], 'r')

            for rec in SeqIO.parse(assembly_in, "fasta"):
                rec.seq = rec.seq.upper()
                contig_ids[rec.id] = contig_id
                contigs[rec.id] = Contig(contig_id, rec, self.telomeres, self.filenames)
                contig_id += 1

            assembly_in.close()

            if len(contigs) == 0:
                log.error(f"Could not load any contigs from {self.assemblyfile}. Is this a valid FASTA file?")
                sys.exit()

        except IOError:
            log.error(f"Can't load assembly from file {self.assemblyfile}!")
            sys.exit()

        log.info(f"Loaded {len(contigs)} contigs from {self.assemblyfile}")
        # Add dictionary of contig IDs to every contig
        for c in contigs:
            contigs[c].contig_ids = contig_ids

        return contigs

    def sample_reads(self):
        if file_exists(self.filenames['sampled_reads']):
            log.info(f"Will use existing {self.filenames['sampled_reads']}")
        else:
            if self.coverage == 0:
                os.symlink(os.path.abspath(self.readfile), self.filenames['sampled_reads'])
                log.info(f"Using all reads as coverage option is 0")
            else:
                log.info(f"Sampling {self.coverage} times coverage of {len(self)/1000000:.1f} Mb assembly from >{self.minreadlength}bp reads in {self.readfile}")

                with gzopen(self.readfile, 'rt') as all_reads, gzopen(self.filenames['sampled_reads'], 'wt', compresslevel=6) as sampled_reads, tqdm(total=self.coverage, leave=False) as pbar:
                    sampled_bases = read_count = times_coverage = 0
                    readset = []

                    for title, sequence, quality in FastqGeneralIterator(all_reads):

                        if len(sequence) < self.minreadlength:
                            continue

                        readset.append(f"@{title}\n{sequence}\n+\n{quality}\n")

                        read_count += 1
                        sampled_bases += len(sequence)
                        new_times_coverage = round(sampled_bases / len(self))

                        if new_times_coverage > times_coverage:
                            print(''.join(readset), file=sampled_reads, end='')
                            readset = []
                            pbar.update()
                            times_coverage = new_times_coverage

                        if times_coverage == self.coverage:
                            break
                    print(''.join(readset), file=sampled_reads, end='')

                log.info(f"Wrote {read_count} reads ({sampled_bases} bases, {times_coverage} times coverage) to {self.filenames['sampled_reads']}")
                if times_coverage < self.coverage:
                    log.warning(f"Only found {times_coverage} times coverage in reads longer than {self.minreadlength}, not {self.coverage} times; consider reducing minimum read length (-l)")


    def make_bam(self, aligntype):
        bam_filename = self.filenames[f'{aligntype}_bam']
        if file_exists(bam_filename, deps=[self.filenames['sampled_reads'], self.filenames['assembly']]):
            log.info(f"Will use existing {bam_filename}")
            return

        if aligntype == 'reads':
           inputfile = self.filenames['sampled_reads']
           align = minimap2[f"-xmap-ont", \
                            "-a", f"-t{self.cores}", self.filenames['assembly'], inputfile]
        elif aligntype == 'contigs':
            inputfile = self.filenames['assembly']
            # -m and -s minimap2 options are only loosely related to alignment length,
            # so use below the true minimum contig alignment threshold; 
            # alignments will be further filtered when loaded into database
            align = minimap2[f"-xasm20", "-D", "-P", \
                             f"-m{self.min_contig_alignment*0.9}", f"-s{self.min_contig_alignment*0.9}", \
                             "-a", f"-t{self.cores}", self.filenames['assembly'], inputfile]
        else:
            log.error(f"Don't know how to align {aligntype}")
            sys.exit()

        # samtools uses cores-1 because -@ specifies additional cores and defaults to 0
        align = align | samtools["sort", f"-@{self.cores-1}", f"-o{bam_filename}"]
        log.info(f"Aligning {aligntype} {inputfile} to assembly")

        try:
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


    def align_to_assembly(self, aligntype):
        self.make_bam(aligntype)
        self.index_bam(aligntype)


    def load_alignments(self):
        alignments = Alignments(self.filenames['alignments'], self.windowsize)
        alignments.load(self.filenames['reads_bam'], self.filenames['contigs_bam'], self.contigs, self.readoutput, self.min_contig_alignment)
        return alignments


    def process_contigs(self):
        log.info(f"Processing {len(self.contigs)} contigs")
        with Pool(self.cores) as p:
            for contig in tapestry_tqdm(p.imap(process_contig, self.contigs.values()), total=len(self.contigs), desc="Processing contigs"):
                self.contigs[contig.name] = contig


    def get_ploidys(self):
        log.info(f"Calculating ploidy estimates")
        fit_ploidy = partial(get_ploidy, ploidy_depths=self.ploidy_depths)
        with Pool(self.cores) as p:
            for contig in tapestry_tqdm(p.imap(fit_ploidy, self.contigs.values()), total=len(self.contigs), desc="Ploidy estimates"):
                self.contigs[contig.name] = contig


    def contig_report(self):
        log.info(f"Generating contig details")

        try:
            with open(f"{self.outdir}/contig_details.tsv", 'wt') as report_file:
                print("Contig\tLength\tGC%\tMedianReadDepth\tStartTelomeres\tEndTelomeres\tStartMeanReadOverhangBases\tEndMeanReadOverhangBases\tUniqueBases\tUnique%\tCategory\tPloidys\tStartConnectors\tEndConnectors", file=report_file)
                for contigname in sorted(self.contigs, key=lambda c: -len(self.contigs[c])):
                    print(self.contigs[contigname].report(self.gc), file=report_file)
        except IOError:
            log.error(f"Could not write contig report")


    def html_report(self):
        log.info(f"Generating Tapestry report")
        env = Environment(loader=FileSystemLoader(report_folder))
        env.globals['include_file'] = include_file
        template = env.get_template('template.html')

        contig_alignments = {c:self.contigs[c].contig_alignments_json() for c in self.contigs}

        contig_coverage = defaultdict(lambda: defaultdict(int))
        for c1 in self.contigs:
            for c2 in self.contigs[c1].contig_coverage:
                contig_coverage[self.contigs[c1].id][self.contigs[c2].id] = self.contigs[c1].contig_coverage[c2]

        read_alignments = {}
        if self.readoutput:
            read_alignments = {c:self.contigs[c].read_alignments for c in self.contigs}

        with open(f"{self.outdir}/{self.basedir}.tapestry_report.html", 'wt') as html_report:
            print(template.render(
                    windowsize = self.windowsize,
                    outdir = self.basedir,
                    assemblyfile = os.path.basename(self.assemblyfile),
                    options = json.dumps(self.options()),
                    contigs = json.dumps([self.contigs[c].json() for c in self.contigs]),
                    read_alignments = json.dumps(read_alignments),
                    contig_alignments = json.dumps(contig_alignments),
                    contig_coverage = json.dumps(contig_coverage),
                    ploidys = json.dumps({c:self.contigs[c].ploidys for c in self.contigs}),
                    median_depth = self.median_depth
                 ),
                 file=html_report)


    def report(self):
        self.contig_report()
        self.html_report()
