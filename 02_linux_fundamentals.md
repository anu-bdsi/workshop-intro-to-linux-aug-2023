# Workshop 2 - Linux basics 

## Learning objectives 

* Understand the standard input and output of Linux
* Be able to use input and output redirection and pipes
* Be able to count, shuf, sort, find, and replace texts in a file 
* Be able to download data from the internet 
* Be able to transfer data between remote and local machine 
* Understand Linux variables and be able to create one 
* Understand and be able to use for loop 

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

## `>` - Output redirection 

We can redirect the output to a file using the `>` symbol. For example:

```sh
head -n 4 > test.txt
```

This redirects the output of the command `head` to a text file, we can use `cat` to view what's in the result file. 

__Make sure to use correct filename, output redirect overwrite everything in an existing file.__

## `>>` - Append output redirection 

If you want to append contents to an existing file, you can use `>>`. Let's try to append something to file `test.txt`:

```sh
ls >> test.txt 
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

We won't use these commands in our course but they are very useful. If you would like to learn more about them, you can read the help manual or watch tutorials online. 

## `grep` - Global regular expression print 

This command is used to search text in a given file or input, it searches for lines containing a match to given text. 

First, let's create a new text file called `linux_unix.txt` and input the following information for us to search later.

```
linux os
unix os
linux is developed from unix
macos is developed from unix
linux and macos are siblings 
```

To search text in a file:

```sh
grep 'pattern' [filename] 
```

__Please search 'linux', 'unix', and 'macos' in file `linux_unix.txt`.__ 

You can use `grep` with pipe to search things in an output of another command:

```sh
ls | grep 'fastq'
```

Some options to search with functionality:

```sh
grep -i # case-insensitive search 
grep -c # print only the count of matching lines 
```

__Please use case-insensitive option to search 'Linux' and 'LinuX'.__

__Please try the `-c` option.__

## `echo` - Display text to the terminal 

The echo command is used to display text or output to the terminal or standard output (stdout). It simply prints the text or variables to the screen. For example, to print a sentence on the screen:

```sh
echo "Hello, world" 
```

You can use either double quotes or single quotes to quote the sentence. 

## `tr` - Translate 

This command allows you to translate or delete characters. It cannot take a file as the argument, so we normally use `tr` with pipes `|` or input redirection `<`.

To translate 6 to 7:

```sh
echo 'file_6' | tr 6 7
```

The `tr` command cannot translate strings, it only works on characters. For example, if I want to replace the name Jiajia with another name:

```sh
echo 'Welcome Jiajia' | tr 'Jiajia' 'David'
```

It gives us:

```
Welcome Dddidd
```

Which is not something we were expecting. That is because it translates each character rather than find the whole thing and replace it. To find and replace strings, we will use another command called `sed`. 

The `tr` command works really well when we want to :

__1. Convert between lower and upper case__

```sh
tr "[:lower:]" "[:upper:]" < linux_unix.txt
```

__2. Convert between space/tab/newline characters etc__

```sh
# convert space characters to newline characters 
tr ' ' '\n' < linux_unix.txt 
```

__3. Delete characters__ 

Here we use the `-d` option of `tr`. 

```sh
echo 'Hello, world' | tr -d 'o'
```

Delete all the space characters:

```sh
tr -d ' ' < linux_unix.txt 
```

Because the `tr` command doesn't work on files directly, you need to use output redirect `>` to save the result. 

## `sed` - Stream editor 

`sed` can perform text manipulation on files. It reads input line by line, performs specified operations on the text, and outputs the modified result. `sed` supports a variety of commands and options to perform a wide range of text transformations.

__Find and replace string:__

The syntax is `sed s/pattern/replacement/g filename`. 

```sh
sed 's/linux/unix/g' linux_unix.txt 
```

For more functionality of `sed`, you can read the help manual. 

# Downloading and Transfering Data 

## `curl` - Download files from the internet

This command allows you to download files from a server using various proticols, including HTTP, HTTPS, FTP, FTPS, SCP, SFTP, and more.

Let's use it to download another fastq file from a FTP address. 

```sh
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz 
```

__Can you explain to me what the `-O` means?__ 

## `scp` and `rsync` - Secure copy 

`scp` allows you to transfer files between a remote machine to a local machine, since we only use our local machine this time we don't need to transfer. But this command is useful when you working on a remote machine. 

`rsync` allows you to transfer directories with everything in it while `scp` only allows file transfer. 

# Linux variables 

In Linux, variables are used to store data, such as text strings, numbers, or arrays, with a name that serves as a reference to that data. The syntax to create a variable is:

```sh
variable_name=value 
```

For example, to create a variable called `name` with the value `John`:

```sh
name="John"
```

To access the value of a variable, we use the `$` sign followed by the variable name:

```sh
echo "My name is $name"
```

Now, try using single quotes for the above command:

```sh
echo 'My name is $name' 
```

What did you get? It printed out `$name` rather than the actual value we refer to. It is because even though both single quotes `''` and double quotes `""` are used to define string literals, they have some differences in behaviours and effects. 

# Single Quotes and Double Quotes

## Single Quotes `''`

Text enclosed in single quotes is treated as a literal string. This means that everything between the single quotes is preserved exactly as written, and no special characters are interpreted. 

For example, if you have a variable `$var` within single quotes, it will be treated as the literal string "$var", not its value.

Single quotes are often used when you want to ensure that the text remains unchanged, without any variable substitution or interpretation of special characters. 

## Double Quotes `""`

Text enclosed in double quotes allows for variable expansion and interpretation of special characters within the string. For example, if you have a variable `$var` within the double quotes, its value will be inserted into the string. 

# Loop statements

A loop is a powerful programming tool that enables you to execute a set of commands repeatedly. There are two types of loops in Linux, for loop and while loop. Here, we will learn `for` loop as it is more frequently used in data analysis. 

## `for` loop

`for` loop allows you to iterate over a list of items and perform actions on each item in the list. The syntax looks like:

```sh
for thing in list_of_things
do 
    command_1 
    command_2 
    ...
done
```

To print each of the elements in the list, we can run: 

```sh
for i in a b c
do 
    echo $i
done
```

or 

```sh
fruits="apple orange banana"

for fruit in $fruits
do
    echo "I like $fruit"
done
```

For loop is very useful when you want to perform the same activity on a list of things. 

# References

* OpenAI - [ChatGPT](https://chat.openai.com/)  
