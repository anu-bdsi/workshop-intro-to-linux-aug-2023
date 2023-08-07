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
cd ~/workshops/variant-calling/
mkdir -p results/sam results/bam results/bcf results/vcf 
```

__Did you notice any difference between the above code and the previous mkdir commands we ran?__

### 3.2. Indexing the reference genome 

The first step of aligning reads to reference genome is to index the genome. Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment. 

Indexing only needs to be run once. Because when you index a reference genome, it will generate a few result files, and the result files store all the indexing information. The generated files are always and will not disappear if you don't delete it. 

The only reason you would want to create a new index is if you are working with a different reference genome or you are using a different alignment tool. 

The command to index a reference genome is: 

```sh
bwa index [genome_file]
```

__Exercise: index the *E. coli* REL606 genome.__ 

What result did you see? Did you get any output files? What are they? 

### 3.3. Aligning reads to reference genome 

In this step, we will use the BWA-MEM algorithm, which is the latest and is generally recommended for high-quality queries as it is faster and more accurate. 

An example command of how to do the alignment is:

```sh
bwa mem [ref_genome] [sample_R1.fastq] [sample_R2.fastq] > [output.sam]
```

__Exercise: please run the alignment on sample SRR2584863.__ 

It takes about 2 minutes to run. After it finishes running, you will get an output SAM file. 

### 3.4. SAM file 

```
@SQ     SN:CP000819.1   LN:4629812
@PG     ID:bwa  PN:bwa  VN:0.7.17-r1188 CL:bwa mem ref-genome/ecoli_rel606.fasta trimmed-fastq/SRR2584863_1.trim.fastq trimmed-fastq/SRR2584863_2.trim.fastq
SRR2584863.3    99      CP000819.1      2067177 60      138M    =       2067947 798     ATCAACAACACGCGTTTATTGGTCTGGCTGATCACCGCCGCCAGGTTGGCGCAGACAAAGGTTTTACCGATTGACGGGCTAACCCCGGTCATCATCAACACATTGTTCTGCGCCTGCATCATCGCGAAGTGCAGGCTG    C@CFFFFFHHHGHIIJJJJJJIIJIIJJJIIJJJJJGGGGGBGHHHFFBEDDDBBDDDDDD:>CBCCC>5;B?CDBBDBBBCCA?BD5<<@AAACCDCDCDDCDDDCDDCDBBDDBD<@>@@CCBDB<@BACCC<38?    NM:i:0  MD:Z:138        MC:Z:28M      AS:i:138        XS:i:0
SRR2584863.3    147     CP000819.1      2067947 60      28M     =       2067177 -798    CAGCGGAAGATCAACAGAATCTTTATCC    IIIIIIIHHHDHHHIDHHF>DDDDA@<@  NM:i:0  MD:Z:28 MC:Z:138M       AS:i:28 XS:i:0
SRR2584863.5    99      CP000819.1      1231976 60      150M    =       1232448 498     TCCGTAAACATTTTAATGTCGTGCTCGAAAGAACGGTGCGTGAACTGCGCGGCGAACCCTGTTTGCAACTGGAAGAGTTTGCACCGACGAAGCAGGAAATTATCTGTTCCCGCTCGTTTGGTGAACGCATCACGGATTATCCGTCGATGC        @@CFBADBFHFHHFIJIGGGGFGHIIIJJJEHIIIGHIJJGIIJGGIIIJJHFDDDBDDDDCDCCCCCDDCCCDDDDACCCDDCDDDDDBD?>CC@DDDDB@ACDCACCCDDBBBB<@D?@BBCCDDDDBBCDB>B>@CDDDDBDDDDC@     NM:i:0   MD:Z:150        MC:Z:26M        AS:i:150        XS:i:0
SRR2584863.5    147     CP000819.1      1232448 60      26M     =       1231976 -498    GCAATTGATGACGGTAATGGATACAC      GHIGHHIIJIIEHHFHHHEDEDF@CC    NM:i:0  MD:Z:26 MC:Z:150M       AS:i:26 XS:i:0
SRR2584863.7    83      CP000819.1      2562378 60      64M     =       2561695 -747    CGGCACGTTAACCTGCTGTTTGATGAGTTTGAACGCTTCTGCCGCGTCCATCGTCGGTACGGAT      ;B?<<<:BABBB?;;3BBCCCB@>B>BDDBDFC@==;CFBFBG<?GC1C8BFF@FFDDDDD@@@        NM:i:0  MD:Z:64 MC:Z:40M        AS:i:64       XS:i:0
SRR2584863.7    163     CP000819.1      2561695 60      40M     =       2562378 747     GGGCCATTCACCACGCAGCCGATAATCGAAACGTCCATCG     ;?<ABDDDHDHFFGGI@GGBBGGEHEGGG9D<GFGCDH)< NM:i:0  MD:Z:40 MC:Z:64M        AS:i:40 XS:i:0
SRR2584863.8    83      CP000819.1      4263438 60      101M    =       4263359 -180    AAGGTGTAGCGCCCCAGGTTACGCAGACGTTCGGCAATCAGGAACAAAATGATCGGCCAGCCCACCAGGAAGCCCAGCGAGTAAATTAAGCCGTCATAGCC CA3?>95555<5DB?4:@<<5@@<?<8?@@@>CCCCCACECECFFDFFFDB:GHHGBHDEEGD9IJIGIGIHEJIGIIIGJJIJJJJJHHHHHFFFFFBC@ NM:i:0  MD:Z:101        MC:Z:150M       AS:i:101        XS:i:0
SRR2584863.8    163     CP000819.1      4263359 60      150M    =       4263438 180     CGCCACCACCACCAGAGAACCACAGGCCGAAAGAATACGAATCGGCCCTTGTTTCAGACGGTAAGAGGCCACATCGGCAAAGGTGTAGCGCCCCAGGTTACGCAGACGTTCGGCAATCAGGAACAAAATGATCGGCCAGCCCACCAGGAA        CCCFFFFFHHGHHJIJIIHJJIIJJJJJJJJJIJICHJGIIDEHGIIIEHCHFDDDCC@CB;AB;@C@BDB@DACDD@BDBDD4::A@C5@BDBBD?CDCC<9@.5A@D??B@@BB>>:AACD?BDDBC@CA>?9@9<?BD@2<<8A8AA     NM:i:0   MD:Z:150        MC:Z:101M       AS:i:150        XS:i:0
```


 



# References
* Data Carpentry - [Data Wrangling and Processing for Genomics](https://datacarpentry.org/wrangling-genomics/index.html)
