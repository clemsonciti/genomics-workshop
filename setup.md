---
layout: page
title: Setup
permalink: /setup/
---

We are going to download some files into our personal "home"
directories on Palmetto from the web.
The data is in the form of a
`.zip` file which can be found at the following URL.

~~~
https://github.com/clemsonciti/genomics-workshop/raw/gh-pages/data/genomics-workshop.zip
~~~

Because we can't run a graphical web browser on the cluster,
we need another way of downloading the `.zip` file.
Fortunately, there is a command-line tool called `wget`
that is meant for exactly this task.

~~~
[nelle@login001 ~]$ wget https://github.com/clemsonciti/genomics-workshop/raw/gh-pages/data/genomics-workshop.zip
~~~
{: .bash}

~~~
[nelle@login001 ~]$ unzip genomics-workshop.zip
Archive:  genomics-workshop.zip
   creating: genomics-workshop/
   creating: genomics-workshop/shell-basics/
   creating: genomics-workshop/shell-basics/data/
   creating: genomics-workshop/shell-basics/data/molecules/
  inflating: genomics-workshop/shell-basics/data/molecules/cubane.pdb
  inflating: genomics-workshop/shell-basics/data/molecules/ethane.pdb
  inflating: genomics-workshop/shell-basics/data/molecules/methane.pdb
  inflating: genomics-workshop/shell-basics/data/molecules/octane.pdb
  inflating: genomics-workshop/shell-basics/data/molecules/pentane.pdb
  inflating: genomics-workshop/shell-basics/data/molecules/propane.pdb
   creating: genomics-workshop/shell-basics/data/pdb/
.
.
.
~~~
{: .output}

Once you have downloaded the `.zip` file,
you may "unzip" it using the `unzip` command:

~~~
[nelle@login001 ~]$ unzip genomics-workshop.zip
~~~
{: .bash}

~~~
[nelle@login001 ~]$ unzip genomics-workshop.zip
Archive:  genomics-workshop.zip
   creating: genomics-workshop/
   creating: genomics-workshop/shell-basics/
   creating: genomics-workshop/shell-basics/data/
   creating: genomics-workshop/shell-basics/data/molecules/
  inflating: genomics-workshop/shell-basics/data/molecules/cubane.pdb
  inflating: genomics-workshop/shell-basics/data/molecules/ethane.pdb
  inflating: genomics-workshop/shell-basics/data/molecules/methane.pdb
  inflating: genomics-workshop/shell-basics/data/molecules/octane.pdb
  inflating: genomics-workshop/shell-basics/data/molecules/pentane.pdb
  inflating: genomics-workshop/shell-basics/data/molecules/propane.pdb
   creating: genomics-workshop/shell-basics/data/pdb/
.
.
.

~~~
{: .output}

Confirm that a new folder `genomics-workshop` has been created
as a result of extracting the `genomics-workshop.zip` file:

~~~
[nelle@login001 ~]$ ls
~~~
{: .bash}

~~~
genomics-workshop  genomics-workshop.zip

~~~
{: .output}

