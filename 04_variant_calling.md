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

__Exercise: please run the alignment on sample `SRR2584863`.__ 

It takes about 2 minutes to run. After it finishes running, you will get an output SAM file. 

### 3.4. SAM file 

SAM stands for Sequence Alignment Map format. It is a TAB-delimited text format consisting of a header section, which is optional, and an alignment section. Header lines start with `@`, while alignment lines do not. Each alignment line has 11 mandatory fields for essential alignment information such as mapping position, and variable number of optional fields for flexible or aligner specific information. 

We can take a look at our SAM file:

__Question: which command do you use to read a file?__ 

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

We can see that there are 2 header lines then following alignment lines. 

__Question: why is there 2 lines have the same sequence name?__ 

### 3.5. BAM file 

BAM is the compressed binary version of SAM. We will convert the SAM file to BAM format using the `samtools` program with the `view` command and specify the output format in BAM using `-b` option.

```sh
samtools view -b [aligned.sam] > [aligned.bam]
```

This takes about 1 minute to run.

__Exercise: compare the size of the SAM and BAM file, how much the file size has decreased?__

### 3.6. Pipe the two steps together 

Normally, to save the disk space we don't keep the SAM files. We can pipe the alignment step and the convert format step together to get the BAM file directly. 

```sh
bwa mem [ref_genome] [sample_R1.fastq] [sample_R2.fastq] | samtools view -b > [aligned.bam]
```

__Exercise: try this method on sample `SRR2584866`.__

### 3.7. Sort BAM file by coordinates 

If you have a look of the previous SAM file, you will find out that the sequences are ordered the same as the fastq file. The reads in fastq files are ordered by the time they go through the sequencing machine. 

In this step, we will sort the BAM file so the sequences are in order of the position they were mapped to the genome from left to right. 

We sort the BAM file using the `sort` command from `samtools`. `-o` tells the command where to write the output. Run this on sample `SRR2584863`. 

```sh
samtools sort -o [sorted.aligned.bam] [aligned.bam]
```

To better understand the sorted result, we can set the output to SAM file as well. Then we can have a look and compare the sorted SAM with the original one. 

```sh
samtools sort -o [sorted.aligned.sam] [aligned.bam]
```

The sorted SAM file should look like this:

```
@HD     VN:1.6  SO:coordinate
@SQ     SN:CP000819.1   LN:4629812
@PG     ID:bwa  PN:bwa  VN:0.7.17-r1188 CL:bwa mem ref-genome/ecoli_rel606.fasta trimmed-fastq/SRR2584863_1.trim.fastq trimmed-fastq/SRR2584863_2.trim.fastq
SRR2584863.72312        2147    CP000819.1      1       60      116H34M =       4629703 4629812 AGCTTTTCATTCTGACTGCAACGGGCAATATGTC   DDDADDDDDEDDEEEDDDDDCDDBDB<?CCDDDD       NM:i:0  MD:Z:34 MC:Z:110M40H    AS:i:34 XS:i:0  SA:Z:CP000819.1,4629697,+,116M34S,60,0;
SRR2584863.236239       99      CP000819.1      1       60      55S95M  =       304     453     GCATGATATTGAAAAAAATATCACCAAATAAAAAACGCCTTAGTAAGTATTTTTCAGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAG        @BCFFFFFHHHHHJJJJJIIJJIJJJJIIJJJJJJIGJIJIIIAHGHCFGIJJJIIIHHHHHGHFFFFFBDECACCDDDB=BBCDDEDDEDCDCDDDB>@@>CDDDDBBDCDDACC@@CCCCC:<A@CCCD@@>C>CC@C>A4<8099?B      NM:i:0  MD:Z:95 MC:Z:150M       AS:i:95 XS:i:0  SA:Z:CP000819.1,4629758,+,55M95S,60,0;
SRR2584863.275624       99      CP000819.1      1       60      14S136M =       145     294     AGTAAGTATTTTTCAGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACTTTAA        @B@FFFBDHHHHHJJIGJJJJJJJJJJJJIJJJJJJIJJJIJIGIGJJJJJJIJFGIHIJJJJJIJHFFFDFEEEEEEEDDDDDDDDDDEDDDDDCDDDDDCBDDBDDDCDDDEEDDDDEDDDEDDDDDDDDDDDDDDDDCDCDDDDEDD      NM:i:0  MD:Z:136        MC:Z:150M       AS:i:136        XS:i:0
SRR2584863.275984       2209    CP000819.1      1       60      87H63M  =       848     896     AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGA       CCCCCCCCCACDDDCCCCCCCCBB>BBBCCCE>>A@A>ACABBC@>>:@>@B@BB<>:@>CCC NM:i:0  MD:Z:63 MC:Z:49M        AS:i:63       XS:i:0  SA:Z:CP000819.1,4629726,+,87M63S,60,0;
SRR2584863.404450       2209    CP000819.1      1       60      91H59M  =       437     481     AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGT   DDDDDDCDDEDDEEEDDDDDDCDDDBDDDDCEEDDCCDDDDDDCDCDDCDDDDBDDCDC     NM:i:0  MD:Z:59 MC:Z:45M        AS:i:59 XS:i:SA:Z:CP000819.1,4629722,+,91M59S,60,0;
SRR2584863.651714       2145    CP000819.1      1       60      95H55M  =       42      175     AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGA       DDDDDDACDEEEEEEDDDDDCDBBD@BDCCDEDDEDDDCD@BDDDDCCDCCB@B? NM:i:0  MD:Z:55 MC:Z:134M       AS:i:55 XS:i:0  SA:Z:CP000819.1,4629718,+,95M55S,60,0;
SRR2584863.702752       163     CP000819.1      1       60      28S122M =       653     802     ATAAAAAACGCCTTAGTAAGTATTTTTCAGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTC        <@BFFFFFHHHHHJIJIJJJHIIJJJJJIJJJJJJHHIIJJJJJJIGIJIJIIJJJJJJJHGHHHHHFDFFFEDECEDDDDBBCACDDDDDEDACDDDDDDDDEDDDCDDDDDDDDBDDDBDDDDDCDEDDDDDCDDDDDCDCDDDDCCD      NM:i:0  MD:Z:122        MC:Z:150M       AS:i:122        XS:i:0
```

This takes about 1 minute to run. 

SAM/BAM files can be sorted in multiple ways, e.g. by location of alignment on the chromosome, by read name, etc. It is important to be aware that different alignment tools will output differently sorted SAM/BAM, and different downstream tools require differently sorted alignment files as input.

For example, the alignment tool BWA-MEM that we used gave us read name sorted file, and the variant calling requires coordinate sorted BAM file as input. 

## 4. Variant calling 

A variant call is a conclusion that there is a nucleotide difference compare to the reference genome at a given position, often referred to as a Single Nucleotide Variant (SNV). The call is usually accompanied by an estimate of variant frequency and some measure of confidence. 

Similar to other steps in this workflow, there are a number of tools available for variant calling. We will use `bcftools`. 

There are a few steps we need to do before the actual call of the variants. 

### 4.1. Calculate the read coverage of positions in the genome

This step is to generate a pileup format representation of sequence read data aligned to a reference genome. A pileup format is a textual representation that shows the aligned reads at each position in the reference genome, along with various summary statistics and variant information. 

```

```


# References
* Data Carpentry - [Data Wrangling and Processing for Genomics](https://datacarpentry.org/wrangling-genomics/index.html)
* The SAM/BAM Format Specification Working Group - [Sequence Alignment Map Format Specification](https://samtools.github.io/hts-specs/SAMv1.pdf) 
* Wikipedia - [Binary Alignment Map](https://en.wikipedia.org/wiki/Binary_Alignment_Map) 
