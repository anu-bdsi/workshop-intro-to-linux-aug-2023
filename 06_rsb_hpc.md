# Workshop 6 - Introduction to the RSB computing cluster and job scheduling system 

## Learning objectives 

* Know the basic information of RSB computer cluster 
* Understand the data storage locations on RSB cluster 
* Understand what is job scheduling system
* Understand the important parts of a SBATCH script 
* Be able to write and run a SBATCH script
* Be able to use parallel processing 

## RSB Computer Cluster 

RSB provides a group of bioinformatics orientated Linux servers for research and teaching. 

The RSB computer cluster consist of 3 servers. The 3 servers work together and are controlled and scheduled by software. 

__dayhoff:__

* 1TB of RAM
* 196 cores
* 100TB of data storage
* Ubuntu 20.04 Linux
* GPU capable (coming soon)

__wright:__

* 1TB of RAM
* 64 cores
* 70TB of data storage
* Ubuntu 20.04 Linux 

__fisher:__

* 1TB of RAM
* 56 cores
* 50TB of data storage 
* Ubuntu 20.04 Linux

## Data Storage Locations

__Home Directory:`/mnt/data/(server)/home/u_id`__ 

You will be automatically in your home directory when you log in to a server. The path to your home directory is `/mnt/data/(server)/home/u_id`. You will have 3 home directories on the 3 servers, and where are you depends on which server you logged into. 

Home directory is where you should install the software to and is where you can store your permanent small data. 

__Groups Directory:`/mnt/data/(server)/home/groups`__

This directory contains directories for every RSB lab. It's for sharing files with lab members. 

__Projects Directory:`/mnt/data/(server)/home/projects`__

This directory contains directories for projects, it is used for sharing and collaborating across groups. 

__Scratch Space:`/mnt/data/(server)/home/scratch`__

The scratch space is for all temporary, working, and non-persistent data. It is allocated 60% of all local storage and is not backed up. Files in the scratch space will be deleted after 130 days. 

A reminder email will be sent out on the first day of each month, containing a list of files that will be deleted at the end of the month. You will need to move important file to somewhere else and delete other files. 

## Job Scheduling System Slurm and PBS Professional

A job scheduling system, often referred to as a workload management or cluster management system, is a software tool designed to efficiently allocate and manage computing resources in a distributed computing environment. 

These systems are commonly used to high-performance computing (HPC) clusters, data centers, and other large-scale computing infrastructures. Their primary purpose is to optimise the utilisation of available resources while ensuring fair access to those resources for multiple users. 

Slurm and PBS Professional are two popular job scheduling systems. They serve similar purpose but different origins and features. 

Slurm is an open-source project that originated from Lawrence Livermore National Library in the United States. It was first developed in the early 2000s and has since gained popularity in the HPC community. PBS Pro has a longer history and was initially developed by NASA in the 1990s. It was later commercialised by several companies, including Altair Engineering, which currently develops and supports PBS Pro. 

Slurm is used on the RSB computer cluster, and NCI uses PBS Pro. 

## Example script for submitting jobs on Slurm 

### `SBATCH` header 

The `SBATCH` header is where we change the allocation settings of a job. We can change parameters like job name, log file paths, time, RAM, CPU numbers etc. 

```sh
#!/bin/bash 
#SBATCH --job-name=job_00
#SBATCH --output=/path/to/output/job_00.out
#SBATCH --error=/path/to/output/job_00.err
#SBATCH --partition=Standard
#SBATCH --time=1:00:00                      # [hh:mm:ss]
#SBATCH --mem=5G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=u_id@anu.edu.au
#SBATCH --mail-type=ALL 
```

* `job-name`: the name of your job. You can name it anything. 
* `output`: path and filename to store the output log file of this job. It contains information that should be printed on the screen.
* `error`: path and filea to store the error log file of this job. It contains all the error messages has in this job. 
* `partition`: is the queue you want to submit your job to. It is used to separating jobs to different queues. We only have one partition on the cluster which is called `Standard`, so all the jobs will be submitted to the same queue. 
* `time`: time limit for this job. Format in `days-hh:mm:ss`. 
* `mem`: RAM size you want to allocate. 
* `nodes`: each node is an independent computer. On our cluster, we have 3 nodes which is dayhoff, wright, and fisher. If you are using parallel processing, you may want to separate your sub-jobs to different nodes to spread the workload. 
* `ntasks`: the number of tasks included in this job. It is used when doing parallel processing with `srun` command. If there is no `srun` command, your entire script would be 1 task. 
* `cpus-per-taks`: number of CPUs you want to allocate for each task. 
* `mail-user`: the email address to receive messages. 
* `mail-type`: types of messages you want to receive. `ALL` for everything. `BEGIN` for job begins execution. `END` for job finishes. `FAIL` for job fails. 

### Activate conda environment in SBATCH script 

The command to activate a conda environment in the SBATCH script is different from the command line. To use a certain conda environment for your script, we need to put the following line of code in the script:

```sh
source /opt/conda/bin/activate /mnt/data/wright/home/[u_id]/.conda/envs/[env-name]
```

### Avoid using `cd` command in the SBATCH script

Because through the job scheduling system, your job can be allocated to any of the nodes. And `cd` command doesn't work across nodes. 

### Write a SBATCH script for the trimming and filtering step

Let's create a script called `run_trim.sh` and input the following code:

```sh
#!/bin/bash 
#SBATCH --job-name=variant_calling 
#SBATCH --output=/mnt/data/wright/home/[u_id]/workshops/variant-calling/variant_calling.out
#SBATCH --error=/mnt/data/wright/home/[u_id]/workshops/variant-calling/variant_calling.err
#SBATCH --partition=Standard
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=email_address
#SBATCH --mail-type=ALL 

source /opt/conda/bin/activate /mnt/data/wright/home/[u_id]/.conda/envs/[env-name]

dir=/mnt/data/wright/home/[u_id]/workshops/variant-calling/raw-fastq
trimmedDir=/mnt/data/wright/home/[u_id]/workshops/variant-calling/trimmed-fastq
NexteraPE=/mnt/data/wright/home/u1133824/.conda/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa

mkdir -p $trimmedDir 

for i in $dir/*_1.fastq.gz
do
base=$(basename $i _1.fastq.gz)

fastq1=$dir/${base}_1.fastq.gz
fastq2=$dir/${base}_2.fastq.gz
trim1=$trimmedDir/${base}_1.trim.fastq.gz
trim2=$trimmedDir/${base}_2.trim.fastq.gz
untrim1=$trimmedDir/${base}_1un.trim.fastq.gz
untrim2=$trimmedDir/${base}_2un.trim.fastq.gz

trimmomatic PE $fastq1 $fastq2 \
                $trim1 $untrim1 \
                $trim2 $untrim2 \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:${NexteraPE}.fa:2:40:15

done
```

To submit the job, we need to run `sbatch run_trim.sh` in the command line. 

## Parallel Processing 

There are several ways to do parallel processing.

### Use multiple CPUs 

Many of the software are built with the ability to use multiple CPUs for running the program. Usually there is an option `-t` or `-threads` for you specify the CPU numbers.

Using Trimmomatic as an example, if we read its help manual there is an option called `-threads` for us to specify the CPU number. We can add this option in and run the job again to compare the time consumed. Let's use 4 CPUs. We need to change the `cpus-per-task` in the SBATCH header accordingly. 

```sh
#!/bin/bash 
#SBATCH --job-name=variant_calling 
#SBATCH --output=/mnt/data/wright/home/[u_id]/workshops/variant-calling/variant_calling.out
#SBATCH --error=/mnt/data/wright/home/[u_id]/workshops/variant-calling/variant_calling.err
#SBATCH --partition=Standard
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=email_address
#SBATCH --mail-type=ALL 

source /opt/conda/bin/activate /mnt/data/wright/home/[u_id]/.conda/envs/[env-name]

dir=/mnt/data/wright/home/[u_id]/workshops/variant-calling/raw-fastq
trimmedDir=/mnt/data/wright/home/[u_id]/workshops/variant-calling/trimmed-fastq
NexteraPE=/mnt/data/wright/home/u1133824/.conda/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa

mkdir -p $trimmedDir 

for i in $dir/*_1.fastq.gz
do
base=$(basename $i _1.fastq.gz)

fastq1=$dir/${base}_1.fastq.gz
fastq2=$dir/${base}_2.fastq.gz
trim1=$trimmedDir/${base}_1.trim.fastq.gz
trim2=$trimmedDir/${base}_2.trim.fastq.gz
untrim1=$trimmedDir/${base}_1un.trim.fastq.gz
untrim2=$trimmedDir/${base}_2un.trim.fastq.gz

trimmomatic PE -threads 4 $fastq1 $fastq2 \
                $trim1 $untrim1 \
                $trim2 $untrim2 \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:${NexteraPE}:2:40:15

done
```

### Using `srun` for parallel running tasks 

`srun` can be used to launch parallel jobs. We need to run trimmomatic on 3 different samples, and in the for loop, we ran them one by one. It is possible to run the 3 trimming process parallel. 

We can do it by adding `srun`, and changing the `--mem`, `--nodes`, and `ntasks` accordingly. 

```sh
#!/bin/bash 
#SBATCH --job-name=variant_calling 
#SBATCH --output=/mnt/data/wright/home/[u_id]/workshops/variant-calling/variant_calling.out
#SBATCH --error=/mnt/data/wright/home/[u_id]/workshops/variant-calling/variant_calling.err
#SBATCH --partition=Standard
#SBATCH --time=1:00:00
#SBATCH --mem=40G
#SBATCH --nodes=3
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=email_address
#SBATCH --mail-type=ALL 

source /opt/conda/bin/activate /mnt/data/wright/home/[u_id]/.conda/envs/[env-name]

dir=/mnt/data/wright/home/[u_id]/workshops/variant-calling/raw-fastq
trimmedDir=/mnt/data/wright/home/[u_id]/workshops/variant-calling/trimmed-fastq
NexteraPE=/mnt/data/wright/home/[u_id]/.conda/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa

mkdir -p $trimmedDir 

for i in $dir/*_1.fastq.gz
do
base=$(basename $i _1.fastq.gz)

fastq1=$dir/${base}_1.fastq.gz
fastq2=$dir/${base}_2.fastq.gz
trim1=$trimmedDir/${base}_1.trim.fastq.gz
trim2=$trimmedDir/${base}_2.trim.fastq.gz
untrim1=$trimmedDir/${base}_1un.trim.fastq.gz
untrim2=$trimmedDir/${base}_2un.trim.fastq.gz

srun --ntasks=1 --mem=10G trimmomatic PE -threads 4 $fastq1 $fastq2 \
                $trim1 $untrim1 \
                $trim2 $untrim2 \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:${NexteraPE}.fa:2:40:15 

done
```

Save and submit the job by running `sbatch run_trim.sh`. 

## SBATCH script for the entire workflow 

```sh
#!/bin/bash 
#SBATCH --job-name=variant_calling 
#SBATCH --output=/mnt/data/wright/home/[u_id]/workshops/variant-calling/variant_calling.out
#SBATCH --error=/mnt/data/wright/home/[u_id]/workshops/variant-calling/variant_calling.err
#SBATCH --partition=Standard
#SBATCH --time=2:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-user=email_address
#SBATCH --mail-type=ALL

source /opt/conda/bin/activate /mnt/data/wright/home/[u_id]/.conda/envs/ecoli-vc

home_dir=/mnt/data/wright/home/[u_id]/workshops/variant-calling

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
NexteraPE=${home_dir}/.conda/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa

bwa index $genome

for i in ${raw_dir}/*_1.fastq.gz
do 
    base=$(basename $i _1.fastq.gz)

    echo "Trimming sample $base"

    fq1=${raw_dir}/${base}_1.fastq.gz
    fq2=${raw_dir}/${base}_2.fastq.gz
    trimmed_fq1=${trimmed_dir}/${base}_1.trim.fastq.gz
    trimmed_fq2=${trimmed_dir}/${base}_2.trim.fastq.gz
    removed_fq1=${trimmed_dir}/${base}_1un.trim.fastq.gz
    removed_fq2=${trimmed_dir}/${base}_2un.trim.fastq.gz

    trimmomatic PE -threads 4 $fq1 $fq2 \
                    $trimmed_fq1 $removed_fq1 \
                    $trimmed_fq2 $removed_fq2 \
                    SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:${NexteraPE}.fa:2:40:15

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

# References 

* Data Carpentry - [Data Wrangling and Processing for Genomics](https://datacarpentry.org/wrangling-genomics/)
* slurm workload manager - [Documentation](https://slurm.schedmd.com/documentation.html) 
* RSB IT Infrastructure Wiki - [RSB bioinformatics server user guide](https://infrawiki.rsb.anu.edu.au/doku.php?id=rsb_it:infrastructure:userdoco:studentservers#running_workloads_and_processes_on_the_servers)
* RONIN - [A simple slurm guide for beginners](https://blog.ronin.cloud/slurm-intro/)
* stack overflow - [parallel but different Slurm job step invocations not working](https://stackoverflow.com/questions/35498763/parallel-but-different-slurm-srun-job-step-invocations-not-working) 