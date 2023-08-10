# Workshop 5 - Result visualisation & Shell script writing 

## Learning objectives 

* 
* 

## 1. Viewing the variant calling results in a genome browser

A genome browser is a graphical interface for displaying information for genomic data. Genome browsers enable users to visualise and browse entire genomes with annotated data, including gene prediction, gene structure, protein, expression, regulation, variation, and comparative analysis.

The [Integrative Genomics Viewer](https://software.broadinstitute.org/software/igv/) (IGV) is a stand-alone browser, which has the advantage of being installed locally and providing fast access. Web-based genome browsers, like [Ensembl](https://www.ensembl.org/index.html) or the [UCSC browser](https://genome.ucsc.edu/), are slower, but provide more functionality.

Before visualising the alignment file, we need to index the sorted BAM file. 

```sh
samtools index [aligned.sorted.bam]
```

__Exercise: please index all 3 samples.__ 

Next, we can start to visualise our data. 

1. Open IGV.
2. Load our reference genome file `ecoli_rel606.fasta` into IGV using the "Load Genomes from File ..." option under the "Genomes" pull-down menu.
3. Load our BAM file `SRR2584866.aligned.sorted.bam` using the "Load from File ..." option under the "File" pull-down menu.
4. Do the same with our VCF file `SRR2584866_final_variants.vcf`.

Your IGV browser should look like:

![igv-result](figures/igv-result.png)

There should be two tracks: one corresponding to our BAM file and the other for our VCF file.

Try zoom in to inspect the variants. 

# References 

* Melbourne Bioinformatics - [Introduction to Genome Browsers](https://www.melbournebioinformatics.org.au/tutorials/tutorials/Genome_browsers/GenomeBrowsers_Intro/)
* Data Carpentry - [Data Wrangling and Processing for Genomics](https://datacarpentry.org/wrangling-genomics/index.html)