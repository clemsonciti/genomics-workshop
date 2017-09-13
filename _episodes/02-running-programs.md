---
title: "Running programs using the shell"
teaching: 15
exercises: 15
objectives:
- "Learn how to run a program using the shell"
- "Learn to redirect standard output to a file, or to another program"
- "Learn to write a loop to run the same program on many input data"
- "Learn to write a shell script"
questions:
keypoints:
---

The shell is more than a way to
create and work with files and folders;
it's also a powerful way to interact with programs.
Most programs that you are familiar with,
for example Microsoft Word or Mozilla Firefox
are *graphical* in nature, i.e., you use
menus, icons, text boxes, and other graphical elements
to control these programs.

Most fields in scientific computing (including computational biology)
are however dominated by command-line programs.

Many programs have **both** a command-line and a graphical interface.

Running and controlling programs via the command-line rather than a GUI
improves the repeatability and reproducibility of your work,
and allows you to automate them.
In this next part of the lesson,
we'll run a simple command-line program (called `wordfreq`).
Given a text file containing several words,
this program counts the frequency of each word,
and prints the result.

Navigate to the directory `genomics-workshop/running-programs/`,
and view the files there:

~~~
$ ls
wordfreq     words-1.txt       words-2.txt       words-3.txt       words-4.txt
~~~

Use the `-l` switch to see a detailed, long list of the files and directories:

~~~
ls -l
~~~

~~~
total 72
-rwxr--r--  1 nelle  cuuser   695 Jul 26 09:14 wordfreq
-rw-r--r--  1 nelle  cuuser  1389 Jul 20 14:40 words-1.txt
-rw-r--r--  1 nelle  cuuser  2573 Jul 20 14:40 words-2.txt
-rw-r--r--  1 nelle  cuuser  8715 Jul 20 14:41 words-3.txt
-rw-r--r--  1 nelle  cuuser  9319 Jul 20 14:41 words-4.txt
~~~

Among other information, such as the
owner of the files (`nelle`), the "user group" that the files belong to (`cuuser`),
the last accessed time, etc.,
the first column of the output shows you
the **permissions** that different entitites have for each file.

A permission of `-rwxr--r--` means:

* The first column `-` means the object is a regular file, for directories,
this column would hold a `d`.
* The next three columns `rwx` tell us what permission the owner of the file has.
In this case, the owner of the file may **read** (`r`), **write** (`w`),
and **execute** (`x`) the file.
While reading/writing applies to any kind of file,
execution only applies to scripts or compiled programs.
Such files are thusly named **executables**.
* The following three columns `r--` tell us what permission the members of
the user group associated with the file have.
In this case, members of the `cuuser` group may **read** (`r`) the file only.
They may not write to it (change its contents).
* The final three columns `---` tell us what permissions everybody else has.
In this case, users that are not part of the `cuuser` group
may not read or write to the file at all.

We see that we have execute `x` permission on the `wordfreq` program.
To run an executable, you can usually just type in
the path to it:

~~~
$ ./wordfreq
~~~
{: .bash}

~~~
usage: wordfreq [-h] [-j N] files [files ...]
wordfreq: error: the following arguments are required: files
~~~
{: .error}

We see that `wordfreq` is missing a required argument
or input `files`.
Let's examine the "help" for `wordfreq` to understand
this better:

~~~
$ ./wordfreq -h
~~~~
{: .bash}

~~~
usage: wordfreq [-h] [-j N] files [files ...]

positional arguments:
  files       Input files

optional arguments:
  -h, --help  show this help message and exit
  -j N        Number of tasks to use
~~~
{: .output}


Let's try specifying the name of an input file:

~~~
$ ./wordfreq words-1.txt
~~~
{: .bash}


~~~
ashwin@laptop ~/w/g/c/g/d/g/running-programs> ./wordfreq words-1.txt 
6 it
9 is
11 a
2 truth
1 universally
1 acknowledged
8 that
4 single
4 man
5 in
2 possession
11 of
~~~
{: .output}

This time, `wordfreq` works as expected,
printing the frequency of each word in the file `words-1.txt`.

So far, all the commands we have been entering either produce no output
(such as `cd`), or print output to the shell window.
It is possible to *redirect* the ouput of a command
to a file using the redirection operator `>`:

~~~
$ ./wordfreq words-1.txt > count-words-1.txt
~~~
{: .bash}

The above command does not print anything.
This is because the output has been *redirected* to the file `count-words-1.txt`.
Let's examine that file to confirm our output is there:

~~~
$ cat count-words-1.txt 
~~~
{: .bash}

~~~
6 it
9 is
11 a
2 truth
1 universally
1 acknowledged
8 that
4 single
4 man
5 in
2 possession
11 of
1 good
~~~
{: .output}

~~~
echo name
echo $name
name=newton
echo $name
~~~

~~~
name=einstein
echo $name
~~~

~~~
for name in einstein newton tesla
do
echo $name
echo name
done
~~~

~~~
for file in *.txt
do
echo $file
done
~~~

~~~
$ ./wordfreq words-1.txt > count-words-1.txt
$ ./wordfreq words-2.txt > count-words-2.txt
$ ./wordfreq words-3.txt > count-words-3.txt
$ ./wordfreq words-4.txt > count-words-4.txt
~~~

~~~
for file in words-*.txt
do
./wordfreq $file > count-$file
done
~~~
