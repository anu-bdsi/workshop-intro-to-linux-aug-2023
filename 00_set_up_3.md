# Set up conda environment on the RSB Wright server

__1. Add conda channels and set priorities__

```sh
conda config --prepend channels conda-forge
conda config --prepend channels bioconda 

conda config --show channels
```

The correct output should look like:

```
channels:
  - conda-forge
  - bioconda
  - defaults
```

__2. Create a conda environment__

```sh
conda create -n ecoli-vc
```

__3. Install packages__

```sh
conda activate ecoli-vc

conda install -c bioconda fastqc
conda install -c bioconda trimmomatic
conda install -c bioconda bwa
conda install -c bioconda samtools=1.9
conda install -c bioconda bcftools=1.9
```

__4. Test if the packages are working__

```sh
fastqc -h
trimmomatic
bwa
samtools
bcftools 
```