# Workshop 4 - Variant calling workflow 

## Learning objectives 

* 
* 

## 3. Alignment to a reference genome 

We mentioned before that we are working with files from a long-term evolution study of an *E. coli* population (designated Ara-3). Now that we have looked at out data to make sure that it is high quality, we can perform variant calling to see how the population changed over time. We care how this population changed relative to the original population *E. coli* strain REL606. Therefore, we will align each of our samples to the *E. coli* REL606 reference genome, and see what differences exist in our reads versus the genome. 

We perform read alignment or mapping to determine where in the genome our reads originated from. There are a number of tools to choose from and, while there is no gold standard, there are some tools that are better suited for particular Next Generation Sequencing (NGS) analyses. 

We will use the [Burrows Wheeler Aligner (BWA)](https://bio-bwa.sourceforge.net/), which is a software package for mapping low-divergent sequences against a large reference genome. It consists of 3 algorithms: BWA-backtrack, BWA-SW, and BWA-MEM. 

BWA-backtrack is designed for Illumina sequence reads up to 100bp. while the rest two for longer sequences ranged from 70bp to 1MBP. BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it is faster and more accurate. BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina reads. 

The alignment process consists of two steps:

1. Indexing the reference genome
2. Aligning the reads to the reference genome 

### 3.1. Downloading the reference genome 

First we need to download the reference genome. 

*Type these step by step, make sure you understand what we are doing each step. Use `tab` key to auto-complete.*

```sh
mkdir -p ~/workshops/variant-calling/ref-genome
cd ~/workshops/variant-calling/ref-genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
mv GCA_000017985.1_ASM1798v1_genomic.fna.gz ecoli_rel606.fasta.gz
gunzip ecoli_rel606.fasta.gz
```

Then, we need to create directories for the results that will be generated as part of this workflow. 

```sh
mkdir -p results/sam results/bam results/bcf results/vcf 
```

__Did you note any difference between the above code and the previous mkdir commands we ran?__

### 3.2. Indexing the reference genome 

The first step of aligning reads to reference genome is to index the genome. Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment. 

Indexing only needs to be run once. Because when you index a reference genome, it will generate a few result files, and the result files store all the indexing information. The generated files are always and will not disappear if you don't delete it. 

The only reason you would want to create a new index is if you are working with a different reference genome or you are using a different alignment tool. 


# References
