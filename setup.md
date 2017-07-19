---
layout: page
title: Setup
permalink: /setup/
---

## Logging-in to the cluster

### Mac OS X and Linux users

Mac OS X or Linux users may open a Terminal, and
type in the following command:

~~~
$ ssh username@login.palmetto.clemson.edu 
~~~
{: .bash}

where `username` is your Clemson user ID.
You will be prompted for both your password and DUO authentication.

### Windows users

Windows users will need to download the MobaXterm package from
[here](<http://mobaxterm.mobatek.net/download.html>).
After downloading and installing MobaXterm:

1.  Launch the MobaXterm program

2.  On the top-left corner of MobaXterm, click the **Session** button. Select
    the SSH setting and confirm that the following settings are set:

    Parameter           |   Value
    --------------------|-------------------------------------
    Remote host         | `login.palmetto.clemson.edu`
    Port                |  22
    X11-Forwarding      | enabled
    Compression         | enabled
    Remote environment  | Interactive shell
    SSH-browser type    | **SCP (enhanced speed)**

3.  Click **OK** and a new session window will be opened, where you will be
    prompted for your Palmetto password and DUO authentication.

After logging in, you should see the following **prompt** on the
last line, indicating that you are logged-in to the `login001`
node of the cluster, and that the shell is waiting for you
type in a command:

~~~
[nelle@login001 ~]$
~~~
{: .bash}

## Downloading workshop data

We are going to download some files into our personal "home"
directories from the web.
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

