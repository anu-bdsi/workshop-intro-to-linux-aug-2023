# Workshop 6 - Introduction to the RSB computing cluster and job scheduling system 

## Learning objectives 

* 
* 

....... 

## Script for submitting jobs in SLURM 

```sh
#!/bin/bash 
#SBATCH --job-name=job_name
#SBATCH --output=/mnt/data/dayhoff/home/u_id/.../job_name.out
#SBATCH --error=/mnt/data/dayhoff/home/u_id/.../job_name.err
#SBATCH --partition=Standard
#SBATCH --time=2:00:00
#SBATCH --mem=5G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=email_address
#SBATCH --mail-type=ALL

HOME_DIR=/mnt/data/dayhoff/home/u_id


```

## Writing a shell script for variant calling 

Create a new file called `variant_calling.sh`, and input the following codes, change the path and emaill address accordingly for the SBATCH headings: 

```sh
#!/bin/bash 
#SBATCH --job-name=variant_calling 
#SBATCH --output=/mnt/data/dayhoff/home/u_id/workshops/variant-calling/variant_calling.out
#SBATCH --error=/mnt/data/dayhoff/home/u_id/workshops/variant-calling/variant_calling.err
#SBATCH --partition=Standard
#SBATCH --time=120
#SBATCH --mem=5G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=email_address
#SBATCH --mail-type=ALL

set -e

source /opt/conda/bin/activate /mnt/data/dayhoff/u_id/.conda/envs/ecoli-vc

home_dir=~/workshops/variant-calling

raw_dir=${home_dir}/raw-fastq
trimmed_dir=${home_dir}/trimmed-fastq
results_dir=${home_dir}/results

mkdir -p ${trimmed_dir} ${results_dir}

sam_dir=${results_dir}/sam
bam_dir=${results_dir}/bam
bcf_dir=${results_dir}/bcf
vcf_dir=${results_dir}/vcf

mkdir -p ${sam_dir} ${bam_dir} ${bcf_dir} ${vcf_dir}

genome=${home_dir}/ref-genome/ecoli_rel606.fasta 

bwa index $genome

for i in ${raw_dir}/*_1.fastq
do 
    base=$(basename $i _1.fastq)

    echo "Trimming sample $base"

    fq1=${raw_dir}/${base}_1.fastq
    fq2=${raw_dir}/${base}_2.fastq
    trimmed_fq1=${trimmed_dir}/${base}_1.trim.fastq
    trimmed_fq2=${trimmed_dir}/${base}_2.trim.fastq
    removed_fq1=${trimmed_dir}/${base}_1un.trim.fastq
    removed_fq2=${trimmed_dir}/${base}_2un.trim.fastq

    trimmomatic PE -threads 4 $fq1 $fq2 \
                    $trimmed_fq1 $removed_fq1 \
                    $trimmed_fq2 $removed_fq2 \
                    SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15

    echo "Aligning sample $base"

    sam=${sam_dir}/${base}.aligned.sam
    bam=${bam_dir}/${base}.aligned.bam
    sorted_bam=${bam_dir}/${base}.aligned.sorted.bam 
    raw_bcf=${bcf_dir}/${base}_raw.bcf 
    variants=${vcf_dir}/${base}_variants.vcf 
    final_variants=${vcf_dir}/${base}_final_variants.vcf 

    bwa mem -t 4 $genome $trimmed_fq1 $trimmed_fq2 | samtools view -S -b > $bam 
    samtools sort -o $sorted_bam $bam 
    samtools index $sorted_bam 

    echo "Calling variants for sample $base"

    bcftools mpileup --threads 4 -O b -o $raw_bcf -f $genome $sorted_bam 
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf 
    vcfutils.pl varFilter $variants > $final_variants 
done
```

Save and exit, now you can run the workflow by:

```sh
sbatch variant_calling.sh 
```

After running the sbatch command, you'll see a job ID popped up on your screen. Wait for a few seconds, and you can run ```squeue``` to check if your job is running. 

__Questions:__ 

* What is ```set -e```?
* Why do we need to add ```/mnt/data/(server)/```?
* What is cluster, node, core, and CPU?

__Exercise:__ 

The samples we just did variant calling on are part of the long-term evolution. The ```SRR2589044``` sample was from generation 5000, ```SRR2584863``` was from generation 15000, and ```SRR2584866``` was from generation 50000. How did the number of mutations change in the sample over time? 

```sh
for file in ~/intro_to_linux/results/vcf/*_final_variants.vcf 
do 
    echo ${file}
    grep -v "#" ${file} | wc -l
done 
```

For ```SRR2589044``` from generation 5000 there were 10 mutations, for ```SRR2584863``` from generation 15000 there were 25 mutations, and ```SRR2584866``` from generation 766 mutations.  

# Use SLURM to run your jobs on RSB IT infrastructure  

Slurm is an open source, fault-tolerant, and highly scalable cluster management and job scheduling system for large and small Linux clusters. It allows users to allocate resources for their jobs rather than taking up whatever resources available. It also let you run jobs in the background rather than looking at the screen all the time. 

## Data storage locations 

There are three nodes on Dayhoff: ```dayhoff```, ```wright```, and ```fisher```. 

If you run ```pwd``` in your home directory, you'll probably see something like ```/home/UID```. But in reality, that's only the path from your home node. For cluster-wide shared data namespace, we have to add ```/mnt/data/(server)``` before the path we get from ```pwd```. 

## A sample sbatch file 

A few things to remember when submitting a SLURM job on the cluster:

* All medium to large processes or workloads need to be run via SLURM, not directly on the command line. If the job was run directly it will be terminated after a short time. 
* With a multinode cluster all the tools you need to use must be installed on all the nodes, and all the file paths need to use the cluster-wide shared data namespace with ```/mnt/data/(server)``` at front. 
* You can limit your job only running on one cluster if you don't want to set up your environment on all nodes. In ```sbatch``` and ```srun```, you can use ```--exclude=fisher, wright``` to remove the nodes you don't want to use. 
* The cluster has only one partition ```standard``` for using right now, but it will have more options in the future for urgent and GPU-intensive workloads. 
* When using multiple cores for a job, first you need to make sure the software you are using supports multi-core processing. Secondly, make sure you specifies the number of cores you want in the codes to run the software as well in the ```sbatch``` file. Sometimes SLURM and software cannot communicate very well so they don't know information from each other. 

A example of what a ```sbatch``` script looks like:

```sh
#!/bin/bash 
#SBATCH --job-name=JobX
#SBATCH --output=/mnt/data/(server)/home/uxxxxx/../%j.%x.out
#SBATCH --error=/mnt/data/(server)/home/uxxxxx/.../%j.%x.err
#SBATCH --partition=Standard
#SBATCH --exclude=wright,fisher 
#SBATCH --time=120:00:00    # 5 days then stop job if not complete
#SBATCH --mem-per-cpu=7G  # 7GB per cpu (rather than per node)
#SBATCH --nodes=2	    # use 2 nodes
#SBATCH --ntasks=Y	    # don't let more than Y tasks run at once
#SBATCH --mem=100G	    # reserve 230GB RAM per node (rather than per cpu)
#SBATCH --cpus-per-task=15  # reserve 15 cpus/threads per task
#SBATCH --ntasks-per-node=Z # only allow z tasks per node
#SBATCH --mail-user uxxxxxxx@anu.edu.au # mail user on job state changes
#SBATCH --mail-type TIME_LIMIT,FAIL		# state changes

srun --ntasks=1 --mem=2G --exclusive my_script_A.sh &
srun --ntasks=1 --mem=2G --exclusive my_script_A.sh &
..
..
srun --ntasks=1 --mem=2G --exclusive my_script_B.sh &
wait
```

It has more options you can use to reserve the resources and manage the job, for more information you can see the [slurm documentation](https://slurm.schedmd.com/overview.html). 

## Parallel Processing 

In the sbatch script above, each srun means a task and it will run parallel with each other. In the variant calling workflow, we used 3 different samples, we can write a script to make the 3 samples run in parallel to save some time. 

Create a file named ```bwa_parallel.sh``` in the directory ```scripts```, and input the following codes, change the path and email address accordingly: 

```sh
#!/bin/bash
#SBATCH --job-name=bwa_p
#SBATCH --output=/mnt/data/dayhoff/home/u_id/intro_to_linux/bwa_p.out
#SBATCH --error=/mnt/data/dayhoff/home/u_id/intro_to_linux/bwa_p.err
#SBATCH --partition=Standard
#SBATCH --exclude=wright,fisher
#SBATCH --time=60
#SBATCH --mem=2G
#SBATCH --nodes=1
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-node=3
#SBATCH --mail-user=email_address
#SBATCH --mail-type=ALL

set -e

indir=~/intro_to_linux/data/trimmed_fastq
outdir=~/intro_to_linux/results
genome=~/intro_to_linux/data/ref_genome/ecoli_rel606.fasta

srun --ntasks=1 --mem=500M bwa mem -t 2 $genome \
            ${indir}/SRR2584863_1.trim.fastq.gz ${indir}/SRR2584863_2.trim.fastq.gz \
            > ${outdir}/sam/SRR2584863.full.aligned.sam &

srun --ntasks=1 --mem=500M bwa mem -t 2 $genome \
            ${indir}/SRR2584866_1.trim.fastq.gz ${indir}/SRR2584866_2.trim.fastq.gz \
            > ${outdir}/sam/SRR2584866.full.aligned.sam &

srun --ntasks=1 --mem=500M bwa mem -t 2 $genome \
            ${indir}/SRR2589044_1.trim.fastq.gz ${indir}/SRR2589044_2.trim.fastq.gz \
            > ${outdir}/sam/SRR2589044.full.aligned.sam &

wait
```

Save ane exit, run ```sbatch bwa_parallel.sh``` to submit the job. You'll see a job ID prompted on your screen. Then, you can run ```squeue``` to check your job status. 

Another useful command to check how much resources have been used for your job is:

```sh
sacct --format=JobID,JobName,State,Start,End,CPUTime,MaxRSS,NodeList,ExitCode --jobs=JOB_ID 
```

## Non-parallel processing 

To compare parallel processing with non-parallel processing, we can write another script to submit the same job and compare the time they used to run the same workload. 

Create a new script called ```bwa_single.sh``` and input the following code, change the path and email address accordingly: 

```sh
#!/bin/bash
#SBATCH --job-name=bwa_s
#SBATCH --output=/mnt/data/dayhoff/home/u_id/intro_to_linux/bwa_s.out
#SBATCH --error=/mnt/data/dayhoff/home/u_id/intro_to_linux/bwa_s.err
#SBATCH --partition=Standard
#SBATCH --exclude=wright,fisher
#SBATCH --time=60
#SBATCH --mem=2G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=email_address
#SBATCH --mail-type=ALL

set -e

indir=~/intro_to_linux/data/trimmed_fastq
genome=~/intro_to_linux/data/ref_genome/ecoli_rel606.fasta

for fq1 in ${indir}/*_1.trim.fastq.gz 
do 
    base=$(basename $fq1 _1.trim.fastq.gz) 

    fq1=~/intro_to_linux/data/trimmed_fastq/${base}_1.trim.fastq.gz
    fq2=~/intro_to_linux/data/trimmed_fastq/${base}_2.trim.fastq.gz
    sam=~/intro_to_linux/results/sam/${base}.aligned.sam

    bwa mem -t 2 $genome $fq1 $fq2 > $sam 
done 

wait
```

Using sbatch to submit the job and compare the result with parallel processing. 

# IGV 

There are a few things we need to use IGV to view the variants:

* The reference genome: ```ecoli_rel606.fasta```.
* The BAM file ```SRR2584866.aligned.sorted.bam``` and the BAM index file ```SRR2584866.aligned.sorted.bam.bai```.
* The VCF file ```SRR2584866_final_variants.vcf```. 

We need to download these files to our laptop and load it through IGV, use ```scp``` to download files from the server. 

After downloading the files, open IGV:

1. Load the reference genome.
2. Load the BAM file. 
3. Load the VCF file. 

Then you can see the variants result. 

# References 

* Data Carpentry - [Data Wrangling and Processing for Genomics](https://datacarpentry.org/wrangling-genomics/)
* slurm workload manager - [Documentation](https://slurm.schedmd.com/documentation.html) 
* RSB IT Infrastructure Wiki - [RSB bioinformatics server user guide](https://infrawiki.rsb.anu.edu.au/doku.php?id=rsb_it:infrastructure:userdoco:studentservers#running_workloads_and_processes_on_the_servers)
* RONIN - [A simple slurm guide for beginners](https://blog.ronin.cloud/slurm-intro/)
* stack overflow - [parallel but different Slurm job step invocations not working](https://stackoverflow.com/questions/35498763/parallel-but-different-slurm-srun-job-step-invocations-not-working) 