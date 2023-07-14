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

__Make sure to use correct filename, output redirect overwrite everything in the file that has the same name.__

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

We won't use these commands in our course but they are very useful. If you would like to learn more about them, you can use the `--help` option to read the manual. 

## `grep` - 