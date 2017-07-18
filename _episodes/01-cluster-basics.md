---
title: "Cluster Basics"
teaching: 20
exercises: 1
start: true
questions:
- "How do I manage files and folders on the cluster?"
- "How do I transfer data to and from the cluster?"
- "How do I submit jobs to the cluster?"
- "Where can I get help?"
keypoints:
- "Keypoint 1"
- "Keypoint 2"
---

## Why use a high-performance computing cluster?

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

After logging in, you should see the following screen:

~~~
 -----------------------------------------------------------------------------
       Welcome to the PALMETTO CLUSTER at CLEMSON UNIVERSITY

   * Email ithelp@clemson.edu with questions or to report problems.

   * Palmetto "office hours" are every Wednesday 8am-11am Barre Hall room B102

   * Quarterly maintenance periods:  late May (followed by Top 500 benchmark),
     early August, late November and early Feb.  Email will be sent
     before each period with details of cluster availability.

        User guide:  http://www.palmetto.clemson.edu/palmetto
   Sample programs:  https://github.com/clemsonciti/palmetto-examples
           Jupyter:  https://www.palmetto.clemson.edu/jupyterhub/docs

   Useful commands:
     module avail             - list available software packages
     qstat -xf jobid          - check status of your job
     qstat -Qf queuename      - check status of a queue
     checkquota               - check your disk quota
     checkqueuecfg            - check general workq max running limits
     cat /etc/hardware-table  - list node hardware: ram,cores,chip,etc.
     qpeek                    - look at a running job's stdout or stderr
     whatsfree                - see what nodes are free right now

   Please do not use /home as your PBS working directory.  Jobs with /home
   as working directory may be killed as performance deteriorates.

   DO NOT RUN JOBS/PROGRAMS/TESTS/PRE-OR-POST PROCESSING ON THE LOGIN NODE.
   They will be terminated without notice. No exceptions.

 -------------- This file is: /etc/motd -------- Last Updated: 25-JUN-2017 ---
~~~
{: .bash}

You will also see the following "prompt", ending with a dollar (`$`) sign:

~~~
[atrikut@login001 ~]$
~~~
{: .bash}

The prompt tells you that you are logged in to the `login001` node of the cluster,
and that you are currently in the "home" directory `/home/nelle`, often shortened
to `~`.

## Organizing code and data on the cluster

The command-line or **shell** is another way of interacting with a computer,
just like the graphical interface (windows, pointer, icons, toolbars and menus) that you may be more familiar with.
You are probably familiar with using a graphical file browser
such as **Windows Explorer** (Windows) or **Finder** (Mac OS X)
for doing the following tasks:

* Navigating between folders
* Creating new files or folders
* Moving, renaming, deleting or copying files and folders

The same tasks can be performed on the command-line,
using a few simple commands:


| Command 	  			| Action
------------------------|-------------------------------------------
| `pwd`		  			| Print working directory
| `cd DIR` 				| Change directory to `DIR`	
| `ls`					| List contents of current directory
| `mkdir DIR`			| Make new directory `DIR`
| `cp SRC DEST`			| Copy file `SRC` into file `DEST`
| `cp -r SRC DEST`		| Copy directory `SRC` into directory `DEST`
| `rm FILE`				| Remove (delete) file `DIR`
| `rm -r DIR`			| Remove (delete) directory `DIR`

We will explore each of these as we go along.

When you login to the cluster,
you begin at the home directory:

~~~
[atrikut@login001 ~]$ pwd
~~~
{: .bash}

~~~
/home/atrikut
~~~
{: .output}

If this is your first time using the cluster,
the output is empty. Let's fix that by creating a directory:

~~~
[atrikut@login001 ~]$ ls
~~~
{: .bash}

~~~

~~~
{: .output}

~~~
[atrikut@login001 ~]$ mkdir genomics-workshop
[atrikut@login001 ~]$ ls
~~~
{: .bash}

~~~
genomics-workshop
~~~
{: .output}

> ## File and folder names
> Notice that we used a hyphen (`-`) between the words
> `genomics` and `workshop` when we created the directory.
> Always avoid using spaces in file and folder names.
> To see why it's a bad idea, try running the following command:
>
> [atrikut@login001 ~]$ mkdir genomics workshop
> 
> Now, type `ls`. What do you see?
>
{: .callout}

## Downloading data from the Web

Before proceeding further,
it is necessary for us to download some files into our home directories
from the web. The data is in the form of a
`.zip` file which can be found at the following URL.

~~~
https://github.com/clemsonciti/genomics-workshop/raw/gh-pages/data/shell-basics.zip
~~~

Because we can't run a graphical web browser on the cluster,
we need another way of downloading the `.zip` file.
Fortunately, there is a command-line tool called `wget`
that is meant for exactly this task.

~~~
[atrikut@login001 ~]$ wget https://github.com/clemsonciti/genomics-workshop/raw/gh-pages/data/shell-basics.zip
~~~
{: .bash}

~~~
--2017-07-18 11:12:53--  https://github.com/clemsonciti/genomics-workshop/raw/gh-pages/data/shell-basics.zip
Resolving github.com (github.com)... 192.30.253.112, 192.30.253.113
Connecting to github.com (github.com)|192.30.253.112|:443... connected.
HTTP request sent, awaiting response... 302 Found
Location: https://raw.githubusercontent.com/clemsonciti/genomics-workshop/gh-pages/data/shell-basics.zip [following]
--2017-07-18 11:12:53--  https://raw.githubusercontent.com/clemsonciti/genomics-workshop/gh-pages/data/shell-basics.zip
Resolving raw.githubusercontent.com (raw.githubusercontent.com)... 151.101.56.133
Connecting to raw.githubusercontent.com (raw.githubusercontent.com)|151.101.56.133|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 58972 (58K) [application/zip]
Saving to: ‘shell-basics.zip’

100%[====================================================================================================================>] 58,972      --.-K/s   in 0.009s

2017-07-18 11:12:53 (6.20 MB/s) - ‘shell-basics.zip’ saved [58972/58972]

~~~
{: .output}

Once you have downloaded the `.zip` file,
you may "unzip" it using the `unzip` command:

~~~
[atrikut@login001 ~]$ unzip shell-basics.zip
~~~
{: .bash}

~~~
Archive:  shell-basics.zip
   creating: shell-basics/
   creating: shell-basics/data/
   creating: shell-basics/data/molecules/
  inflating: shell-basics/data/molecules/cubane.pdb
  inflating: shell-basics/data/molecules/ethane.pdb
  inflating: shell-basics/data/molecules/methane.pdb
  inflating: shell-basics/data/molecules/octane.pdb
  inflating: shell-basics/data/molecules/pentane.pdb
  inflating: shell-basics/data/molecules/propane.pdb
  .
  .
  .
~~~
{: .output}

Confirm that a new folder `shell-basics` has been created
as a result of extracting the `shell-basics.zip` file:

~~~
[atrikut@login001 ~]$ ls
~~~
{: .bash}

~~~
genomics-workshop  shell-basics  shell-basics.zip

~~~
{: .output}

## Filesystem layout

Root directory, home directory
Using relative v/s absolute paths

## Exercise, create a directory structure.

## Types of files

## Using a text editor

## Running scripts

## File and folder permissions

## Getting help
