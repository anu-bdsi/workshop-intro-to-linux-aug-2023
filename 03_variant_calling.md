# Workshop 3 - Variant calling workflow 

## Learning objectives 

* 
* 

# Variant calling 

## What is variant calling?

Variant calling is a bioinformatics process used in genomics to identify genetic variations between the genome of an individual and a reference genome. These variations can include single-nucleotide polymorphisms (SNPs), small insertions and deletions (indels), and larger structural variants. 

In this workshop, our variant calling workflow will only detect SNPs. For finding indels or larger structural variants, the workflow will change a bit. 

![SNPs](figures/snp-quick-ref.jpg)

## The data 

* The data we are going to use is part of a long-term evolution experiment led by [Richard Lenski](https://en.wikipedia.org/wiki/E._coli_long-term_evolution_experiment). 
* The experiment was designed to assess adaptation in E. coli. A population was propagated for more than 40,000 generations in a __glucose-limited__ minimal medium (in most conditions glucose is the best carbon source for E. coli, providing faster growth than other sugars). This medium was supplemented with __citrate__, which E. coli cannot metabolise in the aerobic conditions of the experiment. 
* Sequencing of the populations at regular time points revealed that spontaneous citrate-using variant (__Cit+__) appeared between 31,000 and 31,500 generations, causing an increasing in population size and diversity. 

We will be working with 3 sample events from the Ara-3 strain of this experiment, one from 5,000 generations, one from 15,000 generations, and one from 50,000 generations. The population changed substantially during the course of the experiment, and we will be exploring how with our variant calling workflow. 

We have downloaded the files, and it should be in our `~/workshops/variant-calling/raw-fastq/` directory. Let's go into that directory and have a look.

```sh
cd ~/workshops/variant-calling/raw-fastq/
ls
```

You should see:

![fastq-gz-files](figures/fastq-gz-files.png)

You might notice that our fastq files have the .gz extension behind, it means it's a compressed file and we need to decompress it. The command for decompress a .gz file is `gunzip`. 

To decompress all of our 6 files:

```sh
gunzip *.fastq.gz
```

As we learned the last time, `*` wild card can represent filenames, and we can use it to decompress 6 files in one line of command. It takes several minutes to run. After it finish running, you can use `ls` to check the decompressed files. 

![decompressed](figures/decompressed-files.png)

Then, we can use the `ls -lh` to have a look of how big is our files. 

![file-size](figures/file-size.png)

## Paired-end sequencing 

You might notice that our files contain `_1` and `_2` in the filename, and before that the filenames are the same. It means our sample were sequenced from both ends and it came out as two separate files. In the analysis, we use both of the files. 

## FASTQ file format 

FASTQ format is a text-based format for storing both a biological sequence (usually nucleotide sequence) and its corresponding quality scores. 

It was originally developed to bundle a FASTA formatted sequence and its quality data, but has recently become the de facto standard for storing the output of high-throughput sequencing instruments such as the Illumina Genome Analyzer. 

A FASTQ file has four line-separated fields per sequence:

* Line 1 begins with a `@` character and is followed by a sequence identifier and an optional description.
* Line 2 is the raw sequence letters. 
* Line 3 begins with a `+` character and is optionally followed by the same sequence identifier (and any description) again. 
* Line 4 encodes the quality values for the sequence in line 2, and must contain the same number of symbols as letters in the sequence. 

We can view the first complete read in one of our files by using head command to look at the first 4 lines:

```sh
head -n 4 SRR2584863_1.fastq 
```

And you should see:

![fastq-4-lines](figures/fastq-4-lines.png)

Line 4 shows the quality for each nucleotide in the read, the quality here means the base call accuracy (e.g. 90%). Base calling is the method used to decode the raw signals generated during the sequencing process into the corresponding DNA sequence. 

To make it possible to line up each individual nucleotide with its quality score, the numerical score is converted into a code where it uses letters and special characters to represent the score. There are different encoding format where characters mean different scores, the quality control software usually can automatically detect the encoding format.

# Bioinformatic workflow 

When working with high-throughput sequencing data, the raw reads you get off of the sequencer will need to pass through a number of different tools in order to generate your final desired output. The execution of this set of tools in a specified order is commonly referred to as a workflow or a pipeline. 

There are 9 steps for our variant calling workflow, and we will learn it step by step with hands-on coding exercise. 

![workflow](figures/vc-workflow.png)

## 1. Assessing quality using FastQC 

FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to give a quick impression of whether your data has any problems of which you should be aware before doing any further analysis. 

Rather than looking at quality scores for each individual read, FastQC looks at quality collectively across all reads within a sample. The image below shows one FastQC-generated plot that indicates a very high quality sample:

![fastqc-sample](figures/fastqc-sample.png)

The x-axis displays the base position in the read, and the y-axis shows quality scores. In this example, the sample contains reads that are 40 bp long. This is much shorter than the reads we are working with in our workflow. 

For each position, there is a box-and-whisker plot showing the distribution of quality scores for all reads at that position. The horizontal red line indicates the median quality score and the yellow box shows the 1st to 3rd quartile range. This means that 50% of reads have a quality score that falls within the range of the yellow box at that position. The whiskers show the absolute range, which covers the lowest (0th quartile) to highest (4th quartile) values. 

For each position in this sample, the quality values do not drop much lower that 32. This is a high quality score. The plot background is also colour-coded to identify good (green), acceptable (yellow), and bad (red) quality scores. 

Now let's take a look at a quality plot on the other end of the spectrum. 

![bad-fastqc](figures/bad-fastqc.png)

Here, we see positions within the read in which the boxes span a much wider range. Also, quality scores drop quite low into the "bad" range, particularly on the tail end of the reads. The FastQC tool produces several other diagnostic plots to assess sample quality, in addition to the one plotted above. 

### 1.1. Running FastQC on our data 

First, make sure you have activate your conda environment.

```sh
conda activate ecoli-vc
```

Then, let's validate that FastQC is installed by getting the help document for FastQC.

```sh
fastqc -h 
```

You should see: 

```

            FastQC - A high throughput sequence QC analysis tool

SYNOPSIS

        fastqc seqfile1 seqfile2 .. seqfileN

    fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
           [-c contaminant file] seqfile1 .. seqfileN

DESCRIPTION

    FastQC reads a set of sequence files and produces from each one a quality
    control report consisting of a number of different modules, each one of
    which will help to identify a different potential type of problem in your
    data.

    If no files to process are specified on the command line then the program
    will start as an interactive graphical application.  If files are provided
    on the command line then the program will run with no user interaction
    required.  In this mode it is suitable for inclusion into a standardised
    analysis pipeline.

    The options for the program as as follows:

    -h --help       Print this help file and exit

    -v --version    Print the version of the program and exit

    -o --outdir     Create all output files in the specified output directory.
                    Please note that this directory must exist as the program
                    will not create it.  If this option is not set then the
                    output file for each sequence file is created in the same
                    directory as the sequence file which was processed.

    --casava        Files come from raw casava output. Files in the same sample
                    group (differing only by the group number) will be analysed
                    as a set rather than individually. Sequences with the filter
                    flag set in the header will be excluded from the analysis.
                    Files must have the same names given to them by casava
                    (including being gzipped and ending with .gz) otherwise they
                    won't be grouped together correctly.

    --nano          Files come from nanopore sequences and are in fast5 format. In
                    this mode you can pass in directories to process and the program
                    will take in all fast5 files within those directories and produce
                    a single output file from the sequences found in all files.

    --nofilter      If running with --casava then don't remove read flagged by
                    casava as poor quality when performing the QC analysis.

    --extract       If set then the zipped output file will be uncompressed in
                    the same directory after it has been created. If --delete is
                    also specified then the zip file will be removed after the
                    contents are unzipped.

    -j --java       Provides the full path to the java binary you want to use to
                    launch fastqc. If not supplied then java is assumed to be in
                    your path.

    --noextract     Do not uncompress the output file after creating it.  You
                    should set this option if you do not wish to uncompress
                    the output when running in non-interactive mode.

    --nogroup       Disable grouping of bases for reads >50bp. All reports will
                    show data for every base in the read.  WARNING: Using this
                    option will cause fastqc to crash and burn if you use it on
                    really long reads, and your plots may end up a ridiculous size.
                    You have been warned!

    --min_length    Sets an artificial lower limit on the length of the sequence
                    to be shown in the report.  As long as you set this to a value
                    greater or equal to your longest read length then this will be
                    the sequence length used to create your read groups.  This can
                    be useful for making directly comaparable statistics from
                    datasets with somewhat variable read lengths.

    --dup_length    Sets a length to which the sequences will be truncated when
                    defining them to be duplicates, affecting the duplication and
                    overrepresented sequences plot.  This can be useful if you have
                    long reads with higher levels of miscalls, or contamination with
                    adapter dimers containing UMI sequences.


    -f --format     Bypasses the normal sequence file format detection and
                    forces the program to use the specified format.  Valid
                    formats are bam,sam,bam_mapped,sam_mapped and fastq


    --memory        Sets the base amount of memory, in Megabytes, used to process
                    each file.  Defaults to 512MB.  You may need to increase this if
                    you have a file with very long sequences in it.

    --svg           Save the graphs in the report in SVG format.

    -t --threads    Specifies the number of files which can be processed
                    simultaneously.  Each thread will be allocated 250MB of
                    memory so you shouldn't run more threads than your
                    available memory will cope with, and not more than
                    6 threads on a 32 bit machine

    -c              Specifies a non-default file which contains the list of
    --contaminants  contaminants to screen overrepresented sequences against.
                    The file must contain sets of named contaminants in the
                    form name[tab]sequence.  Lines prefixed with a hash will
                    be ignored.

    -a              Specifies a non-default file which contains the list of
    --adapters      adapter sequences which will be explicity searched against
                    the library. The file must contain sets of named adapters
                    in the form name[tab]sequence.  Lines prefixed with a hash
                    will be ignored.

    -l              Specifies a non-default file which contains a set of criteria
    --limits        which will be used to determine the warn/error limits for the
                    various modules.  This file can also be used to selectively
                    remove some modules from the output all together.  The format
                    needs to mirror the default limits.txt file found in the
                    Configuration folder.

   -k --kmers       Specifies the length of Kmer to look for in the Kmer content
                    module. Specified Kmer length must be between 2 and 10. Default
                    length is 7 if not specified.

   -q --quiet       Suppress all progress messages on stdout and only report errors.

   -d --dir         Selects a directory to be used for temporary files written when
                    generating report images. Defaults to system temp directory if
                    not specified.

BUGS

    Any bugs in fastqc should be reported either to simon.andrews@babraham.ac.uk
    or in www.bioinformatics.babraham.ac.uk/bugzilla/
```

The help document usually gives you guidelines about how to use the software. 

Then, let's change directory to where we store our data files.

```sh
cd ~/workshops/variant-calling/raw-fastq/
```

FastQC can accept multiple files as input, and on both zipped and unzipped files, since we have unzipped all of our files, we can use wild card `*.fastq` to represent our file names. 

```sh
fastqc *.fastq
```

You will see an automatically updating output message telling you the progress of the analysis:

```
Started analysis of SRR2584863_1.fastq
Approx 5% complete for SRR2584863_1.fastq
Approx 10% complete for SRR2584863_1.fastq
Approx 15% complete for SRR2584863_1.fastq
Approx 20% complete for SRR2584863_1.fastq
Approx 25% complete for SRR2584863_1.fastq
Approx 30% complete for SRR2584863_1.fastq
...
```

It should take about 5 minutes to finish running on all 6 files. When the analysis completes, your prompt will return. 

The FastQC program has created several new files within our `raw-fastq` directory. We can use `ls` to have a look. You should see:

![fastqc-results](figures/fastqc-results.png)

For each input FASTQ file, FastQC has created a `.zip` file and a `.html` file. 

Let's take a look at what is the `.html` file about. To open a `.html` file on WSL, you need to open the File Explorer on your Windows system, and on the left pane you will see a penguin icon and a folder for Ubuntu.

![file-explorer](figures/file-explorer.png) 

Then follow the path `/home/jiajia/workshops/variant-calling/raw-fastq/` to go to the `raw-fastq` directory and open the html files. If you are on a Mac computer, you can open your Finder and press `cmd+shift+H` to go to your home directory, then follow the path to find your html files. 

### 1.2. Decoding the other FastQC outputs 







# References

* Wikipedia - [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format) 
* Data Carpentry - [Data Wrangling and Processing for Genomics](https://datacarpentry.org/wrangling-genomics/index.html)
