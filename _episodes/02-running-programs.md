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
usage: wordfreq [-h] fname
wordfreq: error: too few arguments
~~~

~~~
$ ./wordfreq words-1.txt
11 of
11 a
9 is
8 to
8 that
8 the
6 he
6 it
5 his
5 and
5 in
4 be
4 or
4 mr
4 single
~~~

~~~
$ ./wordfreq words-1.txt > count-words-1.txt
~~~

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
