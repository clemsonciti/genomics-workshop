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
Running and controlling programs via the command-line rather than a GUI
improves the repeatability and reproducibility of your work,
and allows you to automate dull tasks.

Navigate to the directory `genomics-workshop/running-programs/`,
and view the files there:

~~~
$ ls
wordfreq     words-1.txt       words-2.txt       words-3.txt       words-4.txt
~~~

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
