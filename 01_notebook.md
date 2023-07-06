# Week 1 - Basic Linux commands 

## Contents today 

* Linux operating systems
* Command Line Interface (CLI)
* The file system 
    * Structure and file path 
    * Moving between directories 
    * Creating, deleting, copying, moving, and renaming files and directories 
    * Editing and viewing files 
* Getting help with a command 
* Using shortcuts 
* Using tmux 

# Linux operating systems 

## What is Linux? 

Linux is a free and open-source operating system that is widely used in various computing devices, from personal computers to servers, mobile devices, embedded systems, and more. It was created by Linux Torvalds in 1991 and is based on the Unix operating system. Popular Linux distributions include Debian, Fedora Linux, and Ubuntu. 

## Linux in biology research 

Linux is extensively used in biology research due to its versatility, robustness, and vast collection of software tools and libraries. 

__Bioinformatics__: Linux offers a wide range of powerful tools and software packages for tasks such as sequence alignment, variant calling, genome assembly, gene expression analysis, next-generation sequencing analysis and more. Popular tools include BLAST, Bowtie, SAMtools, BWA, and Genome Analysis Toolkit (GATK). Linux is also popular in the field of molecular modelling, molecular dynamics simulations, and drug discovery. Software packages like GROMACS, AMBER, AutoDock, and VMD are widely used on Linux for tasks such as protein structure prediction, ligand docking, and virtual screening. 

__High-performance computing (HPC) and cluster computing__: Linux is the dominant operating system in the field of high-performance computing. It is commonly used on supercomputers and computing clusters for running computationally intensive bioinformatics analyses, simulations, and large-scale data processing tasks. 

__Server and web hosting__: Linux-based servers are widely used for hosting biological databases, online bioinformatics tools, and web-based resources for researchers and the scientific community.

# Command Line Interface (CLI)

A command-line interface is a text-based interface used to interact with a computer operating system or software by typing commands into a terminal or command prompt. In a CLI, users communicate with the computer through text commands rather than using graphical user interfaces (GUIs) with windows, menus, and buttons.

In a CLI, users typically enter commands as text strings followed by pressing the Enter/Return key to execute the command. These commands are interpreted by the operating system or software, which then performs the requested actions or provides the desired information. 

__Please open a new terminal__

You may see something like this:

![new-terminal](figures/new-terminal.png)

This is your command line interface, it may look slightly different depending on your device but they share the similar elements.

* __User name__: before the @ sign is your user name. For the example figure, the user name is `jiajia`.
* __Device name__: after the @ sign is your device name. For the example, the device name is `RSB-072750`.
* __Your location__: after the colon and before the $ sign is your location on the device. For the example, the location we are at is `~` which means the home directory. 
* __Command to run__: after the $ sign is where you type the command. 
* __Environment name (optional)__: at the beginning of the command prompt with a bracket, it means the environment you are in. It only shows when you have installed a environment management software on your device, such as Conda or Mamba. For the example, the environment name is `base`. 

Below is also a good example image, you can try the command in there. 

![cli-example](https://res.cloudinary.com/edlitera/image/upload/c_scale/v1564606702/blog/lhp098yhjklkj89.png) 


