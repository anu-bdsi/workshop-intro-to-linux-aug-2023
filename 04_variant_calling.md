# Workshop 4 - Variant calling workflow 

## Learning objectives 

* Understand what is aligning reads to reference genome
* Understand genome indexing 
* Understand SAM/BAM file format 
* Understand the reason of sorting SAM/BAM file by coordinates
* Understand what is read coverage of positions
* Understand VCF/BCF file format 
* Be able to perform each steps of the variant calling workflow 

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

```sh
mkdir -p ~/workshops/variant-calling/ref-genome
cd ~/workshops/variant-calling/ref-genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/
# curl -O for Mac User
mv GCA_000017985.1_ASM1798v1_genomic.fna.gz ecoli_rel606.fasta.gz
gunzip ecoli_rel606.fasta.gz
```

Then, we need to create directories for the results that will be generated as part of this workflow. 

```sh
cd ~/workshops/variant-calling/
mkdir -p results/sam results/bam results/bcf results/vcf 
```

__Q: what does `-p` mean?__

### 3.2. Indexing the reference genome 

The first step of aligning reads to reference genome is to index the genome. Indexing allows the aligner to quickly find potential alignment sites for query sequences in a genome, which saves time during alignment. 

Indexing only needs to be run once. Because when you index a reference genome, it will generate a few result files, and the result files store all the indexing information. The generated files are always there and will not disappear if you don't delete it. 

The only reason you would want to create a new index is if you are working with a different reference genome or you are using a different alignment tool. 

The command to index a reference genome is: 

```sh
bwa index [genome_file]
```

__Please index the *E. coli* REL606 genome.__ 

__Q: Did you get any output files? What are they?__

### 3.3. Aligning reads to reference genome 

In this step, we will use the BWA-MEM algorithm, which is the latest and is generally recommended for high-quality queries as it is faster and more accurate. 

An example command of how to do the alignment is:

```sh
bwa mem [ref_genome] [sample_R1.fastq] [sample_R2.fastq] > [output.sam]
```

__Please run the alignment on sample `SRR2584863`.__ 

It takes about 2 minutes to run. After it finishes running, you will get an output SAM file. 

### 3.4. SAM file 

SAM stands for Sequence Alignment Map format. It is a TAB-delimited text format consisting of a header section, which is optional, and an alignment section. Header lines start with `@`, while alignment lines do not. Each alignment line has 11 mandatory fields for essential alignment information such as mapping position, and variable number of optional fields for flexible or aligner specific information. 

We can take a look at our SAM file: 

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
...
```

We can see that there are 2 header lines then following with alignment lines. 

__Q: why is there 2 lines have the same sequence name?__ 

### 3.5. BAM file 

BAM is the compressed binary version of SAM. We will convert the SAM file to BAM format using the `samtools` program with the `view` command and specify the output format in BAM using `-b` option.

```sh
samtools view -b [aligned.sam] > [aligned.bam]
```

This takes about 1 minute to run.

__Please compare the size of the SAM and BAM file, how much the file size has decreased?__

### 3.6. Pipe the two steps together 

Normally, to save the disk space we don't keep the SAM files. We can pipe the alignment step and the convert format step together to get the BAM file directly. 

```sh
bwa mem [ref_genome] [sample_R1.fastq] [sample_R2.fastq] | samtools view -b > [aligned.bam]
```

__Please try this method on sample `SRR2584866`.__

### 3.7. Sort BAM file by coordinates 

If you have a look of the previous SAM file, you will find out that the sequences are ordered the same as the fastq file. The reads in fastq files are ordered by the time they go through the sequencing machine. 

In this step, we will sort the BAM file so the sequences are in order of the position they were mapped to the genome from left to right. 

We sort the BAM file using the `sort` command from `samtools`. `-o` tells the command where to write the output. Run this on sample `SRR2584863`. 

```sh
samtools sort -o [aligned.sorted.bam] [aligned.bam]
```

This takes about 1 minute to run. 

To better understand the sorted result, we can set the output to SAM file. Then we can have a look and compare the sorted SAM with the original one. 

```sh
samtools sort -o [aligned.sorted.sam] [aligned.bam]
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
...
```

You can see that the lines are ordered by the position in the chromosome rather than the read names. 

SAM/BAM files can be sorted in multiple ways, e.g. by location of alignment on the chromosome, by read name, etc. It is important to be aware that different alignment tools will output differently sorted SAM/BAM, and different downstream tools require differently sorted alignment files as input.

For example, the alignment tool BWA-MEM that we used gave us read name sorted file, and the variant calling requires coordinate sorted BAM file as input. 

## 4. Variant calling 

A variant call is a conclusion that there is a nucleotide difference compare to the reference genome at a given position, often referred to as a Single Nucleotide Variant (SNV). The call is usually accompanied by an estimate of variant frequency and some measure of confidence. 

Similar to other steps in this workflow, there are a number of tools available for variant calling. We will use `bcftools`. 

There are a few steps we need to do before the actual call of the variants. 

### 4.1. Calculate the read coverage of positions in the genome

This step is to generate a pileup format representation of sequence read data aligned to the reference genome. A pileup format is a textual representation that shows the aligned reads at each position in the reference genome, along with various summary statistics and variant information. 

We will use the `mpileup` command from `bcftools`. The flag `-O v` tells bcftools to generate a VCF format output, `-o` specifies where to write the output file, and `-f` tells the path to the reference genome. 

The Variant Call Format (VCF) is a text file format (most likely stored in a compressed manner). It contains meta-information lines, a header line, and then data lines each containing information about a position in the genome. The format also has the ability to contain genotype information on samples for each position. 

```sh 
bcftools mpileup -O v -o [sample_raw.vcf] -f [reg_genome] [aligned.sorted.bam]
```

It takes about 10 minutes to run. Your result VCF file should look like this:

```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##bcftoolsVersion=1.9+htslib-1.9
##bcftoolsCommand=mpileup -O v -o vcf/SRR2584863_raw.vcf -f ../ref-genome/ecoli_rel606.fasta bam/SRR2584863.sorted.aligned.bam
##reference=file://../ref-genome/ecoli_rel606.fasta
##contig=<ID=CP000819.1,length=4629812>
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##INFO=<ID=I16,Number=16,Type=Float,Description="Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h">
##INFO=<ID=QS,Number=R,Type=Float,Description="Auxiliary tag used for calling">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  bam/SRR2584863.sorted.aligned.bam
CP000819.1      1       .       A       <*>     0       .       DP=27;I16=1,18,0,0,676,24074,0,0,1140,68400,0,0,0,0,0,0;QS=1,0;MQSB=1;MQ0F=0  PL      0,57,229
CP000819.1      2       .       G       <*>     0       .       DP=27;I16=2,18,0,0,745,28387,0,0,1200,72000,0,0,26,68,0,0;QS=1,0;MQSB=1;MQ0F=0        PL      0,60,255
CP000819.1      3       .       C       <*>     0       .       DP=27;I16=2,18,0,0,736,27594,0,0,1200,72000,0,0,46,140,0,0;QS=1,0;MQSB=1;MQ0F=0       PL      0,60,255
CP000819.1      4       .       T       <*>     0       .       DP=27;I16=2,18,0,0,741,27889,0,0,1200,72000,0,0,66,252,0,0;QS=1,0;MQSB=1;MQ0F=0       PL      0,60,255
CP000819.1      5       .       T       <*>     0       .       DP=27;I16=2,18,0,0,754,28670,0,0,1200,72000,0,0,86,404,0,0;QS=1,0;MQSB=1;MQ0F=0       PL      0,60,255
CP000819.1      6       .       T       <*>     0       .       DP=27;I16=2,18,0,0,778,30440,0,0,1200,72000,0,0,106,596,0,0;QS=1,0;MQSB=1;MQ0F=0      PL      0,60,255
CP000819.1      7       .       T       <*>     0       .       DP=27;I16=2,18,0,0,763,29341,0,0,1200,72000,0,0,126,828,0,0;QS=1,0;MQSB=1;MQ0F=0      PL      0,60,255
CP000819.1      8       .       C       <*>     0       .       DP=27;I16=2,18,0,0,786,30962,0,0,1200,72000,0,0,146,1100,0,0;QS=1,0;MQSB=1;MQ0F=0     PL      0,60,255
...
```

Where it stores the information start from the first base of the reference genome. 

* `CHROM` - chromosome
* `POS` - position: the reference position, with the 1st base having position 1. 
* `ID` - identifier
* `REF` - reference base
* `ALT` - alternate base 
* `QUAL` - quality score 
* `FILTER` - filter status: PASS if this position has passed all filters (a call is made at this position).
* `INFO` - additional information 

Like the alignment steps, we don't normally need to read the coverage information so it would be good we can also output this information to a binary file to reduce the storage space used. 

__Please read the documentation of `bcftools mpileup` to see if you can make the output into a binary format? If can please do so.__ 

### 4.2. Detect the single nucleotide variants (SNVs)

The above step we only calculated the read coverage but didn't do anything about it, so it contains everything on every position of the genome. As we can see from the VCF file, most of the positions don't have an alternate base. 

In this step, we will use `bcftools call` to identify SNVs. We need to specify ploidy with the flag `--ploidy`, `-m` allows for multiallelic and rare-variant calling, `-v` to output variant sites only, and `-o` specifies where to write the output file.

```sh
bcftools call --ploidy 1 -m -v -o [sample_variants.vcf] [sample_raw.bcf]
```

The variants VCF file should look like this:

```
...
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  bam/SRR2584863.sorted.aligned.bam
CP000819.1      9972    .       T       G       225     .       DP=81;VDB=0.387587;SGB=-0.693147;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=0,0,31,48;MQ=60  GT:PL   1:255,0
CP000819.1      263235  .       G       T       228     .       DP=65;VDB=0.00108727;SGB=-0.693146;RPB=0.0433688;MQB=3.80355e-08;MQSB=0.368447;BQB=0.00167695;MQ0F=0.246154;AC=1;AN=1;DP4=10,6,15,27;MQ=28    GT:PL   1:255,0
CP000819.1      281923  .       G       T       225     .       DP=69;VDB=0.0125476;SGB=-0.693147;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=0,0,40,24;MQ=60 GT:PL   1:255,0
CP000819.1      377000  .       T       G       228     .       DP=60;VDB=0.575301;SGB=-0.693146;RPB=0.467702;MQB=1.73175e-07;MQSB=0.239114;BQB=0.319368;MQ0F=0.233333;AC=1;AN=1;DP4=12,2,22,21;MQ=30 GT:PL   1:255,0
CP000819.1      433359  .       CTTTTTTT        CTTTTTTTT       116     .       INDEL;IDV=73;IMF=0.948052;DP=77;VDB=0.0026231;SGB=-0.693147;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=1,3,25,48;MQ=60       GT:PL   1:143,0
CP000819.1      473901  .       CCGC    CCGCGC  228     .       INDEL;IDV=68;IMF=0.944444;DP=72;VDB=0.802322;SGB=-0.693147;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=2,1,21,48;MQ=60        GT:PL   1:255,0
CP000819.1      648692  .       C       T       225     .       DP=84;VDB=0.151313;SGB=-0.693147;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=0,0,58,21;MQ=60  GT:PL   1:255,0
CP000819.1      1331794 .       C       A       225     .       DP=54;VDB=0.330798;SGB=-0.693147;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=0,0,27,25;MQ=60  GT:PL   1:255,0
CP000819.1      1733343 .       G       A       225     .       DP=95;VDB=0.996142;SGB=-0.693147;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=0,0,47,44;MQ=60  GT:PL   1:255,0
CP000819.1      2103887 .       ACAGCCAGCCAGCCAGCCAGCCAGCCAGCCAG        ACAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAGCCAG     7.30276  .       INDEL;IDV=12;IMF=0.444444;DP=27;VDB=0.969575;SGB=-0.686358;MQSB=0.951472;MQ0F=0;AC=1;AN=1;DP4=9,4,4,10;MQ=58    GT:PL1:255,221
CP000819.1      2333538 .       AT      ATT     228     .       INDEL;IDV=78;IMF=0.962963;DP=81;VDB=0.283362;SGB=-0.693147;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=0,9,40,32;MQ=60        GT:PL   1:255,0
...
```

Now, there are only variant sites left. 

### 4.3. Filter and report the SNV variants

The above step gave us the result of all variant sites in our sample, but not all the sites have high quality score. In this step, we will filter out low quality variants using `vcfutils.pl varFilter`.

The `vcfutils.pl varFilter` filters out variants that do not meet minimum quality default criteria, which can be changed through its options. 

```sh
vcfutils.pl varFilter [sample_variants.vcf] > [sample_final_variants.vcf]
```

Let's take a look of our final result, it should only contain the filtered variants.

__Exercise: take a guess of how many variants left?__

__Exercise: how can you count it using the vcf file?__

Here, we have finished all the steps for variant calling. Next, we are going to use a genome browser to look at our data. 

# Homework

Run the variant calling workflow on the other 2 samples. 

# References
* Data Carpentry - [Data Wrangling and Processing for Genomics](https://datacarpentry.org/wrangling-genomics/index.html)
* samtools - [Sequence Alignment Map Format Specification](https://samtools.github.io/hts-specs/SAMv1.pdf) 
* Wikipedia - [Binary Alignment Map](https://en.wikipedia.org/wiki/Binary_Alignment_Map) 
* samtools - [The Variant Call Format (VCF) Version 4.2 Specification](http://samtools.github.io/hts-specs/VCFv4.2.pdf)