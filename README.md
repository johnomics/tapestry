# Tapestry

Tapestry is a tool to validate and edit small eukaryotic genome assemblies using long sequence reads. It is designed to help identify complete chromosomes, symbionts, haplotypes, complex features and errors in close-to-complete genome assemblies.

Tapestry takes as input a genome assembly, a set of long reads, and a telomere sequence, and produces a HTML report summarising the assembly. The report can be used to sort, filter and describe contigs based on the summary information. A new set of contigs can be exported from the report and used to filter the original assembly.

The report is intended to be shared with collaborators to explain the genome assembly, or perhaps included as supplemental material for a genome paper.

Tapestry is designed for use with small eukaryotic genome assemblies, perhaps less than 50 Mb and less than 100 contigs. It will run on larger genomes (particularly with the -n option, see below), but it may be slow.

## Getting started

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
- networkx
- numpy
- pandas
- plumbum
- pysam
- scikit-learn >= 0.20
- sqlalchemy
- tqdm


### Basic usage

For a genome with assembly `assembly.fasta`, read set `reads.fastq.gz` and telomere sequence `TGATGA`, run Tapestry like this:

```
weave -a assembly.fasta -r reads.fastq.gz -t TGATGA -o assembly -c 4
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
                        (default 10000)
  -t TELOMERE [TELOMERE ...], --telomere TELOMERE [TELOMERE ...]
                        telomere sequence to search for
  -w WINDOWSIZE, --windowsize WINDOWSIZE
                        window size for ploidy calculations (default 10000)
  -n, --noreadoutput    do not output read alignments in report (default
                        False)
  -o OUTPUT, --output OUTPUT
                        directory to write output, default weave_output
  -c CORES, --cores CORES
                        number of parallel cores to use (default 1)
  -v, --version         report version number and exit
```

`weave` will subsample 50 times coverage of the genome from the FASTQ of long reads, using reads at least 10,000 bases long. It will use this read set to analyse the genome. The coverage can be altered with `-d` and the minimum read length with `-l`.

More than one telomere sequence can be provided, eg `-t TGATGA CTTATT`. (If you don't have a telomere sequence for your genome, try spot-checking the first and last few kilobases of the contigs in your FASTA file; the telomeres are usually quite obvious if the assembly contains complete chromosomes).

`weave` estimates ploidy counts for windows of 10,000 bases by default. This can be adjusted with `-w`.

By default, Tapestry reports contain a plot of read alignments for each contig. However, if your genome assembly is large, the read alignments can make the report large and difficult to render. Turn off the read alignment output with `-n`.

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



## Tapestry report

A Tapestry report is a HTML file that is intended to work in any modern browser (please raise an issue if it does not). The report contains a summary of the assembly, a plot of each contig with associated summary statistics in a neighbouring table, and a plot of read alignments for a selected contig.

### Contig plot and table
The contig plot shows each contig with regions shaded green. The depth of green reflects an estimate of the ploidy for each region (see read alignment plot for key). If telomeres have been found, they will appear as red circles at the end of each contig. The opacity of the circles reflects the number of telomeres found.

The contig plot can be sorted by any column in the contig table, including contig length, GC content and read depth. The table includes a Cluster column, which shows Tapestry's attempt to cluster contigs into groups based on read and contig alignments. These clusters can be edited by clicking on the cluster text.

Contigs can be removed from the assembly by unchecking the checkboxes in the Keep column of the contig table. As contigs are rejected, the assembly table will be updated to show the size and N50s of the remaining contigs and the rejected contigs. By sorting the table by the Keep column first, all of the rejected contigs will be sorted to the bottom of the plot.

Notes can be added for each contig, so contigs can be annotated as complete, haplotypes, symbionts and so on. This column could also be used to impose a particular sort order on the assembly.

Above the contig table, there are buttons to turn off certain columns in the table (for example, if the report is too wide for your screen), to sort the table by multiple columns, and to export the table to a CSV file.

The radio buttons to the left of the table select one contig to show in the read alignments plot (see below).

### Contig alignments
Clicking on a contig name to the left of the contig plot will show contig alignments for that contig. The currently selected contig is shown with a black dashed line; contigs with alignments to the selected contig will be connected with purple lines. The opacity of a contig's purple line reflects the percentage of the selected contig that aligns to that contig. Strong lines indicate the selected contig may be mostly or completely contained in other contigs.

Individual contig alignments are shown as light blue blocks underneath the contig plots. Hovering over an alignment will show the destination of this alignment. Overlapping alignments will increase the blue shading; deep blue regions have alignments to many other contigs.

### Read alignments

The read alignments plot shows the sampled read alignments for one selected contig. Alignments are shown in blue, with left and right clipped regions shown in several colours. If a read has no other alignment to the left or right of this one, the clipped region is shown in beige. If a neighbouring part of the read does align elsewhere, the distance to that alignment is shown in pink (alignment to the same contig) or purple (alignment to a different contig). Hovering over alignments will show the contigs for any neighbouring alignments, or '-' if there are no neighbours. Alignments are shaded by mapping quality; light alignments are low quality.

The green lines in the background show ploidy read depths, estimated by assuming that the median read depth for all windows across the genome represents diploid read depth.

The vertical black dashed lines show the start and end of the contig. Completed chromosomes should show few overlaps beyond these lines.

This plot is intended to give an overview of contig quality, showing whether contigs are complete chromosomes (by read alignments ending at contig ends) or whether they have misassemblies (by breaks in coverage or many alignments to other contigs). It does not replace a proper read alignment viewer; BAM files are available in the output folder for loading into [IGV](https://software.broadinstitute.org/software/igv/) or a similar viewer.

### Export and import

The export button to the top right of the contig table will save the current table to a CSV file called `<output>_filtered.csv`. The CSV file can be imported back into the report using the file chooser. You can therefore save your editing progress at any time and come back to it later. The CSV file can be used to filter the original assembly FASTA with `clean`, or it can be delivered with the report file to show a cleaned up assembly in the report.

### contig_details.tsv

The `contig_details.tsv` file output by Tapestry contains the basic information from the report and some additional details. It reports the following fields for each contig:

- Cluster, Contig, Length, GC%, MedianReadDepth - as per the HTML report. See Connectors below for a description of the clustering process.
- StartTelomeres, EndTelomeres - the number of telomeres found in the 1 kilobase sequence at the start or end of the contig (used to show red circles in the report contig plot)
- StartMeanReadOverhangBases, EndMeanReadOverhangBases - the number of bases overhanging the ends of the contig, based on clipping of reads. For example, if a 10kb read has 2-10kb aligning to 1-9kb of a contig, with 1kb clipped, this read has a 1kb overhang. Short overhangs at the ends of a contig may indicate the contig is a complete chromosome.
- UniqueBases, Unique% - the number of bases with no alignments to other contigs. Low uniqueness often reflects repetitive content or misassemblies; however, even complete chromosomes can have low uniqueness if the genome is highly repetitive.
- Category - a crude attempt to classify contigs. This is a code with three characters; the first is `N` for nuclear sequence, `-` otherwise; the second is an indicator of completeness, `C` if both ends are complete, `L` if the left end is, `R` if the right end is, or `-` otherwise; the third is the majority ploidy number. So `NC2` might be a complete diploid nuclear chromosome; `N-1` might be a haplotype. A contig is labelled nuclear (`N`) if its GC content is within 2% of the mean GC content for the entire assembly. A contig end is labelled complete (`C`,`L`,`R`) if it has at least one telomere sequence and its mean read overhang length is less than 250bp.
- Ploidys - a comma-separated list of percentages of the contig called for each ploidy. For example, `1:0.20,2:0.70,6:0.10` indicates that 20% of the contig is haploid, 70% is diploid and 10% has six copies (repetitive). Ploidys are estimated from read alignments using a Bayesian Gaussian Mixture model. The ploidy with the largest percentage is reported in the Category.
- StartConnectors, EndConnectors - if a contig has an incomplete end, Tapestry looks for read and contig alignments from the 10kb at the contig end that connect to other contigs. If it finds connecting regions, it counts the reads aligned to the original region in total (O), to the destination region in total (D), and the reads connecting the two (C). It then calculates percentages for C/O and C/D, the percentage of reads in the original and destination reads that connect. If C/O is > 65% and C/D > 10%, the connector is reported here, with the two percentages. Contigs with connectors are clustered together.


## Contact

If you have any problems with or comments about Tapestry, please raise an issue or contact [John Davey](mailto:john.davey@york.ac.uk). Thank you!


## Acknowledgements

Tapestry was written by John Davey ([email](mailto:john.davey@york.ac.uk), [Twitter](http://twitter.com/johnomics)) in the [Genomics and Bioinformatics](http://york.ac.uk/biology/technology-facility/genomics/) lab at the University of York, supported by pump priming funding from the [Department of Biology](https://www.york.ac.uk/biology/) to [Seth Davis](https://www.york.ac.uk/biology/research/plant-biology/seth-davis/) and [Jeremy Mottram](https://www.york.ac.uk/biology/research/infection-immunity/mottram/).  
