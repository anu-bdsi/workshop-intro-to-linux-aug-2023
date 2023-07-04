# For Windows users

## 1. Install WSL 

Windows Subsystem for Linux (WSL) is a compatibility layer introduced by Microsoft that allows users to run a Linux environment directly on Windows operating systems. 

Before installing, make sure you have administrator rights for your laptop. Then, you can follow the first 4 steps on this [page](https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-11-with-gui-support#1-overview) to install Ubuntu WSL on your laptop. You don't need to do the 5th step "Install and use GUI package" and 6th step. 

## 2. Install Anaconda in WSL 

Anaconda is an open-source package management system and environment management system commonly used in data science, scientific computing, and machine learning projects. It allows you to install, update, and manage software packages and dependencies for your projects. It supports packages written in various programming languages, with a focus on Python packages. 

First, let's start a Ubuntu terminal. 

Then, run the following codes in your terminal:

```sh
# change directory to home 
cd ~

# wget is a command to download files from the internet
wget https://repo.anaconda.com/archive/Anaconda3-2023.03-1-Linux-x86_64.sh 

# run the installer to install 
bash Anaconda3-2023.03-1-Linux-x86_64.sh 
```

You will be prompted with some questions in the installation process, just press `enter` or type `yes` then `enter` to continue. When you see `--MORE--` at the bottom of the terminal, press `space` key to move to the next page of the document. When you are asked to decide the installation location, just press `enter` to use the default path. 

After you finish the installation, use `exit` command to exit the terminal and restart a new one.

In the new terminal, run the following command to check if Conda has been successfully installed:

```sh
conda --version
```

You should see `conda 23.3.1` printed on screen. 

## 3. Install necessary packages using conda



# For MacOS users 

## 1. Install Anaconda

Anaconda is an open-source package management system and environment management system commonly used in data science, scientific computing, and machine learning projects. It allows you to install, update, and manage software packages and dependencies for your projects. It supports packages written in various programming languages, with a focus on Python packages. 

* If you have an ANU managed device

You can download Anaconda from the university package manager. 

* If you have a personal laptop

Go to this [website](https://www.anaconda.com/download#downloads) and pick the suitable installer for your device, you can pick either graphical or command line installer. I recommend people without command line experience to use the graphical installer. 

## 2. Install necessary packages using conda

