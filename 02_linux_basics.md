# Workshop 2 - Linux basics 

## Contents today 

* 
* 

# Standard input and output in Linux 

In Linux, when you enter a command as an input, you receive an output. It's the basic workflow of Linux. 

By default: 
* The standard input (stdin) device is the keyboard.
* The standard output (stdout) device is the screen. 

# Input and output redirection 

## `<` - Input redirection 

The input can be redirected from a file or another command's output. To redirect input from a file, we use the `<` symbol followed by the filename. For example:

```sh
cat < sample_1.fastq 
```

works the same as `cat sample_1.fastq` which prints the file contents on the screen. 

And there is a command called `mail` in Linux which can send emails, we can redirect a file as the email content to this command rather than type everything in the command line. 

```sh 
mail -s "Subject" email@address < email.txt 
```

__Please create the `email.txt` file by yourself.__

## `>` - Output redirection 

We can redirect the output to a file using the `>` symbol. For example:

```sh
ls -alh /etc > things_in_etc.txt 
```

This redirects the output of the command `ls` to a text file, we can use less to view what's in the result file. 

__Make sure to use correct filename, output redirect overwrite everything in a existing file.__

## `>>` - Append output redirection 

If you want to append contents to an existing file, you can use `>>`. Let's try to append something to file `things_in_etc.txt`:

```sh
ls >> things_in_etc.txt 
```

## `|` - Pipes

Pipe is another form of redirection, it can redirect the output of the first command to the input of the second command. For example: 

```sh
# list only the first 5 files in /etc
ls /etc | head -n 5
```

We can pipe as many commands as we like. For example:

```sh
# list only the fifth file in /etc
ls /etc | head -n 5 | tail -n 1
```

# Text manipulation in Linux 

## `wc` - Word count 

This command is used to count the number of lines, words, and characters in a file. 

```sh
wc sample_1.fastq
```

The output should be:

```
73428   165213 11164775 sample_1.fastq 
```

It is useful when we trying to get an idea of how many data we have inside a file. For example, fastq files are formatted as 4 lines representing one sequence. So, we can know that there are 73428/4=18357 sequences in the sample fastq file. 

There are options to only display one data:
* ```wc -l``` only shows the number of lines. 
* ```wc -w``` only shows the number of words. 
* ```wc -c``` only shows the number of characters. 

## `shuf`, `sort`, and `uniq`

`shuf` allows you to shuffle a file by lines. 

`sort` allows you to sort a file by lines.

`uniq` allows you to remove duplicates of presorted files. 

We won't use these commands in our course but they are very useful. If you would like to learn more about them, you can use the `--help` option to read the manual or watch tutorials online. 

## `grep` - Global regular expression print 

This command is used to search text in a given file or input, it searches for lines containing a match to given text. 

To search text in a file:

```sh
grep 'conf' things_in_etc.txt 
```

You can use `grep` with pipe to search things in an output of another command:

```sh
ls | grep 'fastq'
```

Some options to search with functionality:

```sh
grep -i # case-insensitive search 
grep -v # invert matching, select non-matching lines 
grep -c # print only the count of matching lines 
grep -n # prefix the output with line number 
```

## `tr` - Translate 

This command allows you to translate or delete characters. It cannot take a file as input, so we normally use `tr` with pipes `|`. For example:

To translate 6 to 7:

```sh
echo 'file_6.txt' | tr 6 7
```

The `tr` command cannot translate strings, it only works on characters. For example, if I want to replace the name with another name:

```sh
echo 'Welcome Jiajia!' | tr 'Jiajia' 'David'
```

It gives us:

```
Welcome Dddidd! 
```

Which is not something we were expecting. That is because it translates each character rather than find the whole thing and replace it. To find and replace strings, we will use another command called `sed`. 

The `tr` command works really well when we want to :

__1. Convert between lower and upper case__

```sh
tr [:lower:] [:upper:] < things_in_etc.txt 
```

__2. Convert between space/tab/newline characters etc__

```sh
# convert space characters to newline characters 
tr ' ' '\n' < things_in_etc.txt 
```

__3. Delete characters__ 

Here we use the `-d` option of `tr`. 

```sh
echo 'Hello, world!' | tr -d 'o'
```

Delete all the space characters:

```sh
tr -d ' ' < things_in_etc.txt 
```

Because the `tr` command doesn't work on files directly, you need to use output redirect `>` to save the result. 

## `sed` - Stream editor 

`sed` can perform text manipulation on files. It reads input line by line, performs specified operations on the text, and outputs the modified result. `sed` supports a variety of commands and options to perform a wide range of text transformations.

__Find and replace string:__

The syntax is `sed s/pattern/replacement/g filename`. To find replace the string `root` with `admin` on file `things_in_etc.txt`:

```sh
sed 's/root/admin/g' things_in_etc.txt 
```

For more functionality of `sed`, you can use `--help` to read the manual. 

# Downloading and transfering data 

## `wget` - Download files from the internet

This command allows you to download files from http/https/ftp addresses, we will use it to download the data for the variant calling workflow. 

First, let's create a directory to store the data files.

```sh
mkdir -p ~/workshops/variant-calling/raw-fastq/
cd ~/workshops/variant-calling/raw-fastq/ 
```

The `-p` option of `mkdir` allows to create parent directories as needed. In other words, if any directories in the path don't exist, the command will create them. 

```sh
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_2.fastq.gz
```

## `scp` - Secure copy 

This command allows you to transfer files from a remote machine to a local machine, since we are using our local machine we don't need file trnasfer. But this command is useful when you working on a remote machine. 
