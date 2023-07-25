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

* 
* 


# The workflow 

There are 9 steps in this workflow and we will learn it step by step with hands-on coding exercise. 

![workflow](figures/vc-workflow.png)



