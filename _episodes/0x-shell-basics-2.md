
The shell is more than a way to create and look at files and folders,
it's also a way for you to run programs and view their output.
In fact, so far, we have been doing exactly this.
When you type in `pwd`, you are really running a program called `pwd`.
The program produces some output (in this case, the current working directory),
which you view on the screen.
You can use the `which` command to see where the program is located:

~~~
$ which pwd
/usr/bin/pwd
~~~
{: .bash}

Programs can be **executables**, i.e.,
binary instructions that are ready to be executed by the computer.
Such executables are often created by compiling code written in a programming language
such as C, C++, Java or Fortran.
A special program called a *compiler* is used for this purpose.

Programs can also be stored as **scripts**.
These text files that are
read by a program called an *interpreter*,
which converts the script into instructions for the computer.

Inside the `programs` folder,
we have two scripts, a Python script called `hello.py` and a Bash script called `hello.sh`.
We also some C source-code (`hello.c`), which we will compile into an executable.

~~~
[atrikut@login001 shell-basics]$ cd programs
[atrikut@login001 programs]$ ls
~~~
{: .bash}

~~~
hello.c  hello.py hello.sh
~~~
{: .output}

### Running a Python script

~~~
[atrikut@login001 programs]$ python hello.py
~~~
{: .bash}

~~~
Hello from Python!
~~~
{: .output}


### Running a bash script

~~~
[atrikut@login001 programs]$ bash hello.sh
~~~
{: .bash}

~~~
Hello from Bash!
~~~
{: .output}

### Running a C program

~~~
[atrikut@login001 programs]$ gcc hello.c -o hello.out
[atrikut@login001 programs]$ ./hello.out
~~~
{: .bash}

~~~
Hello from C!
~~~
{: .output}
