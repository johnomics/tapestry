# Tapestry

Tapestry is a tool to validate and edit small eukaryotic genome assemblies using long sequence reads. It is designed to help identify complete chromosomes, symbionts, haplotypes, complex features and errors in close-to-complete genome assemblies.

Tapestry takes as input a genome assembly, a set of long reads, and a telomere sequence, and produces a HTML report summarising the assembly. The report can be used to sort, filter and describe contigs based on the summary information. A new set of contigs can be exported from the report and used to filter the original assembly.

The report is intended to be shared with collaborators to explain the genome assembly, or perhaps included as supplemental material for a genome paper.

Tapestry is designed for manual editing of small eukaryotic genome assemblies, perhaps less than 50 Mb and less than 500 contigs. It will run on larger assemblies, but it will be slow, and some features will be disabled.

## Getting Started

### Installation
Tapestry is in the [Bioconda](https://bioconda.github.io) repository. If you are set up to use Bioconda (see [Install conda](https://bioconda.github.io/#install-conda) and [Set up channels](https://bioconda.github.io/#set-up-channels)), install the latest release with `conda`:

```
conda install tapestry
```

If you want to install the latest commit, make sure the requirements below are satisfied (or a previous release is installed with conda), and then pull from github and install:
```
git clone https://github.com/johnomics/tapestry
cd tapestry
python setup.py install
```

### Requirements

Tapestry requires the following packages and tools (which will be installed by `conda` automatically):

- Linux or macOS
- Python 3.6 or later

- [minimap2](https://github.com/lh3/minimap2)
- [samtools](http://www.htslib.org/download/)

Python packages:
- biopython
- intervaltree
- jinja2
- numpy
- pandas
- plumbum
- pysam
- sqlalchemy
- tqdm


### Basic usage

For a genome with assembly `assembly.fasta`, read set `reads.fastq.gz` and telomere sequence `TTAGGG`, run Tapestry like this:

```
weave -a assembly.fasta -r reads.fastq.gz -t TTAGGG -o assembly -c 4
```

`-o` is an output directory of your choice, `-c` is the number of cores to use (more is recommended).

`weave` will produce a directory of output files, including a `<output>.tapestry_report.html` file - `assembly.tapestry_report.html` in this case. Open this file in a browser of your choice to view and edit the assembly.

After filtering and sorting your assembly, you can export a CSV file called `<output>_filtered.csv` from the report, which will describe the new assembly. You can then produce a cleaned assembly FASTA file with `clean`:

```
clean -a assembly.fasta -c assembly_filtered.csv
```


## weave

`weave` generates an assembly report and several intermediate output files that can be used for further analyses. It has the following options:
```
  -h, --help            show this help message and exit
  -a ASSEMBLY, --assembly ASSEMBLY
                        filename of assembly in FASTA format (required)
  -r READS, --reads READS
                        filename of long reads in FASTQ format (required; must
                        be gzipped)
  -d DEPTH, --depth DEPTH
                        genome coverage to subsample from FASTQ file (default
                        50)
  -l LENGTH, --length LENGTH
                        minimum read length to retain when subsampling
                        (default 10000 bp)
  -t TELOMERE [TELOMERE ...], --telomere TELOMERE [TELOMERE ...]
                        telomere sequence to search for
  -w WINDOWSIZE, --windowsize WINDOWSIZE
                        window size for ploidy calculations (default ~1/30th
                        of contig N50 length, minimum 10000 bp)
  -f, --forcereadoutput
                        output read alignments whatever the assembly size
                        (default, only output read alignments for <50 Mb
                        assemblies)
  -m MINCONTIGALIGNMENT, --mincontigalignment MINCONTIGALIGNMENT
                        minimum length of contig alignment to keep (default
                        2000 bp)
  -o OUTPUT, --output OUTPUT
                        directory to write output, default weave_output
  -c CORES, --cores CORES
                        number of parallel cores to use (default 1)
  -v, --version         report version number and exit
```

`weave` will subsample 50 times coverage of the genome from the FASTQ of long reads, using reads at least 10,000 bases long. It will use this read set to analyse the genome. The coverage can be altered with `-d` and the minimum read length with `-l`.

More than one telomere sequence can be provided, eg `-t TTAGGG CTTATT`. (If you don't have a telomere sequence for your genome, try spot-checking the first and last few kilobases of the contigs in your FASTA file; the telomeres are usually quite obvious if the assembly contains complete chromosomes).

`weave` estimates ploidy counts for windows across the genome. Window size is roughly 1/30th of contig N50 size by default. This can be adjusted with `-w`.

By default, Tapestry reports contain a plot of read alignments for each contig. However, if your genome assembly is large, the read alignments will make the report large and difficult to render, so the read alignment plot is not included in the report for assemblies over 50 Mb in size. The output of read alignments can be forced with `-f`, but reducing read depth with `-d` is strongly recommended.

`weave` generates the following output files in the output directory provided with `-o`:

- `assembly.fasta` - a copy of the assembly
- `reads.fastq.gz` - the subsampled read set
- `contigs_assembly.bam`, `contigs_assembly.bam.bai` - BAM file and index of the assembly mapped to itself with minimap2
- `reads_assembly.bam`, `reads_assembly.bam,bai` - BAM file and index of the subsampled reads mapped to the assembly with minimap2
- `alignments.db` - an SQLite database of read and contig alignments
- `<output>.tapestry_report.html` - the main Tapestry report, to be opened in a web browser
- `contig_details.tsv` - a TSV file containing information about each contig, including some gory details not included in the main report (see below)

## clean

`clean` will filter your original genome assembly based on a CSV exported from the Tapestry report:

```
  -h, --help            show this help message and exit
  -a ASSEMBLY, --assembly ASSEMBLY
                        filename of assembly in FASTA format
  -c CSV, --csv CSV     Tapestry CSV output
  -o OUTPUT, --output OUTPUT
                        filename of output contigs, default
                        filtered_assembly.fasta
```



## Tapestry Report

A Tapestry report is a HTML file that is intended to work in any modern web browser (please raise an issue if it does not). The report contains a summary of the assembly, a plot of each contig with associated summary statistics in a neighbouring table, and a plot of read alignments for a selected contig.

There is an example report, [`c_merolae_simulated.tapestry_report.html`](example/c_merolae_simulated.tapestry_report.html) in the example folder of this repository. The contig annotation file [`c_merolae_simulated.filtered.csv`](example/c_merolae_simulated.filtered.csv) can be loaded into the report using the `Choose File` button at the top of the report. The report is on a [canu 1.9](https://github.com/marbl/canu) assembly of a [badread](https://github.com/rrwick/Badread)-simulated read set of the [C. merolae](https://plants.ensembl.org/Cyanidioschyzon_merolae/Info/Index) genome (20 nuclear chromosomes, 16.5 Mb long, plus mitochondrial and chloroplast genomes), which was modified to include 1 translocation between chromosomes 17 and 19, one triploid chromosome, 20 copy number variants, ~0.1% SNPs and ~0.01% indels using [simuG](https://github.com/yjx1217/simuG). The following should be read while viewing the example report.

The report begins with an assembly table showing basic statistics for the whole assembly, and for sets of selected and rejected contigs, depending on which contigs have been selected in the contigs table (see below).

### Contig plot and table
The contig plot shows each contig with regions shaded green. The depth of green reflects an estimate of the ploidy for each region (see read alignment plot for key). In the example report, tig00000083 is estimated to be diploid, whereas tig00000092 is the triploid chromosome.

If telomeres have been found, they will appear as red vertical lines at the ends of each contig. The opacity of the circles reflects the number of telomeres found; any count above 20 will have maximum opacity. tig00000083 has telomeres at both ends, but tig00003756 only has a telomere at the start.

The contig plot can be sorted by any column in the contig table, including contig length, GC content and read depth. These fields can often be used to identify accessory genomes or to filter out bad quality contigs. In the example report, the chloroplast (tig00000091) and mitochondrial (tig00000004) genomes can be clearly seen by GC content and read depth, and a range of junk contigs (generated due to badread simulating junk nanopore reads) can be rejected due to extreme GC content and low read depth.

The contig table also includes Group and Notes columns. The Group column can be used to cluster together contigs that may be of one type, or that may come from one chromosome or organism; the Notes column can be used to annotate particular contigs. In the example annotation CSV, the groups are the original C. merolae chromosomes (plus a few additional groups for accessory genomes and extraneous contigs), and the notes indicate if the contigs are complete chromosomes, fragmented, haplotypes, or something else. Both of these columns can be sorted, so they could be used to specify a particular order on the contigs that can't be achieved by sorting other columns.

Contigs can be removed from the assembly by unchecking the checkboxes in the Keep column of the contig table. As contigs are rejected, the assembly table will be updated to show the size and N50s of the remaining contigs and the rejected contigs. By sorting the table by the Keep column first, all of the rejected contigs will be sorted to the bottom of the plot. In the example annotation, unique material has been retained, and haplotypes and junk removed.

Above the contig table, there are buttons to turn off certain columns in the table (for example, if the report is too wide for your screen), to sort the table by multiple columns, and to export the table to a CSV file.

The radio buttons to the left of the table select one contig to show in the read alignments plot (see below).

### Contig alignments
Clicking on a contig name to the left of the contig plot will show contig alignments for that contig. The currently selected contig is shown with a black dashed line; contigs with alignments to the selected contig will be connected with purple lines. The opacity of a contig's purple line reflects the percentage of the selected contig that aligns to that contig. Strong lines indicate the selected contig may be mostly or completely contained in other contigs.

Individual contig alignments are shown as light blue blocks underneath the contig plots. Hovering over an alignment will show the destination of this alignment. Overlapping alignments will increase the blue shading; deep blue regions have alignments to many other contigs.

For example, clicking on tig00000099's name in the contig plot shows that there are short alignments from this contig to several other contigs, but long alignments covering the whole contig to tig00003760 (the purple line connecting these two contigs is strong). Also, ploidy is estimated to be haploid or none for the aligned regions of both contigs. So tig00000099 is likely to be a haplotype of tig0000003760 caused by a copy number variation.

### Read alignments

The read alignments plot shows the sampled read alignments for one contig, selected by radio button to the left of the contig table. Alignments are shown in blue, with left and right clipped regions shown in several colours. If a read has no other alignment to the left or right of this one, the clipped region is shown in beige. If a neighbouring part of the read does align elsewhere, the distance to that alignment is shown in grey (alignment to the same contig) or purple (alignment to a different contig). Hovering over alignments will show the contigs for any neighbouring alignments, or '-' if there are no neighbours. Alignments are shaded by mapping quality; light alignments are low quality.

The green lines in the background show ploidy read depths, estimated by assuming that the median read depth for all windows across the genome represents diploid read depth.

The vertical black dashed lines show the start and end of the contig. Completed chromosomes should show few overlaps beyond these lines.

This plot is intended to give an overview of contig quality, showing whether contigs are complete chromosomes (by read alignments ending at contig ends) or whether they have misassemblies (by breaks in coverage or many alignments to other contigs). It does not replace a proper read alignment viewer; BAM files are available in the output folder for loading into [IGV](https://software.broadinstitute.org/software/igv/) or a similar viewer.

For example, tig000003757 and tig00000094 are two halves of chromosome 14 which have failed to assemble complete. tig00003757 has a telomere at its end but not its start. Clicking on the radio button next to tig00003757 will show the read alignment plot for this contig at the bottom of the report. There are many reads that end exactly at the right end of the contig, where the telomere is found, indicating this is a chromosome end. But there are also reads that have overhangs beyond the start of the contig (purple lines to the left of the plot). Hovering over these reads shows that the overhanging regions align to tig00000094. The contig alignments for these contigs also show an overlap, suggesting that these contigs may come from the same chromosome, or have a common repeat.

### Export and import

The export button to the top right of the contig table will save the current table to a CSV file called `<output>.filtered.csv`. The CSV file can be imported back into the report using the file chooser at the top of the report. You can therefore save your editing progress at any time and come back to it later. The CSV file can be used to filter the original assembly FASTA with `clean`, or it can be delivered with the report file to show a cleaned up assembly in the report.

### contig_details.tsv

The `contig_details.tsv` file output by Tapestry contains the basic information from the report and some additional details. It reports the following fields for each contig:

- Contig, Length, GC%, MedianReadDepth - as per the HTML report.
- StartTelomeres, EndTelomeres - the number of telomeres found in the 1 kilobase sequence at the start or end of the contig (used to show red vertical lines in the report contig plot)
- StartMeanReadOverhangBases, EndMeanReadOverhangBases - the number of bases overhanging the ends of the contig, based on clipping of reads. For example, if a 10kb read has 2-10kb aligning to 1-9kb of a contig, with 1kb clipped, this read has a 1kb overhang. Short overhangs at the ends of a contig may indicate the contig is a complete chromosome.
- UniqueBases, Unique% - the number of bases with no alignments to other contigs. Low uniqueness often reflects repetitive content or misassemblies; however, even complete chromosomes can have low uniqueness if the genome is highly repetitive.
- Ploidys - a comma-separated list of percentages of the contig called for each ploidy. For example, `1:0.20,2:0.70,6:0.10` indicates that 20% of the contig is haploid, 70% is diploid and 10% has six copies (repetitive). Ploidys are estimated from read alignments.


## Contact

If you have any problems with or comments about Tapestry, please raise an issue or contact [John Davey](mailto:john.davey@york.ac.uk). Thank you!


## Acknowledgements

Tapestry was written by John Davey ([email](mailto:john.davey@york.ac.uk), [Twitter](http://twitter.com/johnomics)) in the [Genomics and Bioinformatics](http://york.ac.uk/biology/technology-facility/genomics/) lab at the University of York, supported by pump priming funding from the [Department of Biology](https://www.york.ac.uk/biology/) to [Seth Davis](https://www.york.ac.uk/biology/research/plant-biology/seth-davis/) and [Jeremy Mottram](https://www.york.ac.uk/biology/research/infection-immunity/mottram/).  
