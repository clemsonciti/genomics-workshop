---
title: "Interacting with the Palmetto Cluster"
teaching: 15
exercises: 15
objectives:
- "Learn the structure of a typical HPC cluster"
- "Learn about the different filesystems on the cluster"
- "Learn to use built-in software"
- "Learn to submit jobs to the Palmetto cluster"
questions:
- "How is using HPC different from using my laptop or workstation?"
- "What does the cluster \"look\" like?"
- "Where can I store data on the cluster?"
- "What software is available on the cluster?"
- "How do I reserve hardware on the cluster?"
keypoints:
---

Now that you are familiar with using the command-line
to work with files and folders, run programs,
and with writing simple shell scripts;
you can learn a little about using an HPC cluster like Palmetto.

## Differences between HPC and personal laptop/workstation

Using an HPC cluster is in many ways different from using
computing devices you are used to using
such as your laptop, workstation, or mobile phone.

1. Unlike your personal laptop or workstation,
where you may be running an operating system
such as Microsoft Windows, Mac OS X or Ubuntu,
which offer both a graphical interface *and* a command-line interface,
an HPC cluster gives you only a command-line interface.

1. Unlike your personal laptop or workstation,
which you are generally sitting in front of when working with,
you will use the Palmetto cluster **remotely**
over a network.

1. Unlike your personal laptop or workstation,
the cluster is not a single machine,
but rather several small machines
connected together in a network.
Each one of these machines is not very much more
powerful than your own laptop or workstation.
In fact, the vast majority of the world's "supercomputers"
are not single, large, powerful machines, but rather many small machines.

1. Unlike your personal laptop or workstation,
which has the CPU, memory (RAM), and storage
all together in one "box",
the cluster has centrally located storage disks
that all computers are connected to.
This way, when you create a folder called `work`
in your `/home/username` directory on one machine,
the same folder `work` is visible and accessible by all other machines.

1. Unlike your personal laptop or workstation,
on which there is likely only one person running programs
at any given time,
the cluster is a **shared** resource on which
several people are running programs on the cluster at any given time.
The next difference is an important consequence of this.

1. Unlike your personal laptop or workstation,
on which you probably have *administrative priveleges*,
i.e., you can create or delete any files or folders
or install any software you like,
you have very limited permissions on the cluster.
There are only a handful of folders you can
write and store data on (e.g., your personal home directory `/home/username`).
Additionally, you can only install software packages into your own home directory,
which sometimes means that you need to compile the package yourself.
If you have never compiled software before,
this is an additional skill that you will need to learn -- it's not as complicated
as it sounds!

1. Unlike you personal laptop or workstation,
on which you run programs *interactively*,
i.e., you ask your computer to run a program
and it immediately starts to work on it.
On the cluster, you will primarily be computing in **batches**.
Each "batch" or **job** is a set of tasks or commands
that you "submit" to the cluster.
The tasks will run whenever the compute resources required
to run them become available.
Depending on the compute resources required,
your tasks may start immediately,
or may wait in "queue" to run for upto hours or even days.

## Structure of the Palmetto cluster

Here's a picture of the cluster:

<img src="{{site.baseurl}}/fig/palmetto-front-view.png" style="width:500px">

Each one of the "cabinets" in the picture above is a
"rack", containing several computers or **nodes**:

<img src="{{site.baseurl}}/fig/palmetto-nodes-closeup.png" style="width:500px">

The nodes are all connected together in a network.
The vast majority of nodes on the cluster are
**compute** nodes (labeled `node0001`, `node002`, etc.)
which run users' computing tasks and
do most of the heavy-lifting on the cluster.
Currently, the cluster has over 2000 compute nodes.
Users may **not** login to any of the compute nodes directly.

A few nodes are service nodes which
are responsible for other things than intensive computations;
for example:

* the "Login" node (labeled `login001`) is the entry-point to the cluster,
that all users generally first log-in to.
From this node, users may organize their projects
and submit tasks to the cluster.

* the "Data Transfer" node (labeled `xfer01`) is used
to schedule data transfers between Palmetto and other machines.
Users may login to this node when moving data in and out
of the cluster, but they can not submit tasks from this node.

* the "Scheduler" node (labeled `pbs02`) "schedules" tasks
submitted by users to the compute nodes.
Users may **not** log-in to this node.

<img src="{{site.baseurl}}/fig/palmetto-structure.png" style="width:500px">

Finally, while all nodes have their own local storage,
there are also **shared** storage systems
that all nodes are connected to.
The most prominent of these are the storage for users' home directories.
Other shared filesystems are the "scratch" directories used for temporary data storage.
Access to files and folders on a shared storage system
happens over the network,
so it is generally slower than local storage,
but it is extremely convenient for data that needs to be visible
to all nodes,
or for data that needs to be accessed after a job completes.

## Submitting your first job to the cluster

The login node (`login001` or `login.palmetto.clemson.edu`)
is not meant for running any computationally intensive programs
(i.e., programs that take a long time to run,
or use a lot of CPU cores or memory).
We will simulate such a program with the following command:

~~~
$ sleep 100
~~~

(Hit `Ctrl+C` to interrupt the program and regain control of the shell).

> ## Making a sleeper script
>
> Create a script called `sleeper.sh`
> which contains the `sleep` command above.
> Add the following line after the `sleep` command:
>
> echo "Good morning!"
>
> On the login node, run the script using:
>
> ~~~
> $ sh sleeper.sh
> ~~~
>
> And kill it immediately using `Ctrl+C`!
{: .challenge}

We now have an example of a long-running task.
To "submit" this task to the cluster,
we will prepare a "batch script".
A batch script is a regular shell script with a few additions
(for example, a description of the amount of computing resources required to run the commands in the script).
This batch script is then "submitted" to
the **scheduler**.
The scheduler manages all submitted jobs
and allocates compute resources to them as they become available.
Until the resources you request are available,
the scheduler may "queue" your job.

We can convert our shell script `sleeper.sh` to a
batch script by adding the following two lines to the
**top** of the script:

~~~
#PBS -N sleep
#PBS -l select=1:ncpus=1:mem=1gb,walltime=00:05:00
~~~

Normally, lines beginning with a `#` are ignored by the shell
(often useful to add *comments* to scripts).
Lines beginning with `#PBS` are instructions or *directives* for the scheduler.
The line `#PBS -N sleep` tells the scheduler that we'd like to name our job
`sleep` (you can choose any name you like).
The line `#PBS -l select=1:ncpus=1:mem=1gb,walltime=00:01:00`
tells the scheduler to allocate 1 "chunk" of hardware for this job
-- each "chunk" containing 1 CPU core and 1gb of RAM -- for a maximum
walltime of five minutes.

The reason the directives begin with the characters `#PBS`
is that the scheduler software used on the cluster
is called PBS (**P**ortable **B**atch **S**ystem).
Different clusters may run different schedulers,
so the specific commands to submit and control jobs on those clusters
may differ slightly.

> ## Submitting the `sleeper.sh` batch script
>
> After adding the lines above to your script
> `sleeper.sh`, you can submit it to the scheduler
> using the `qsub` command:
>
> ~~~
> $ qsub sleeper.sh
> ~~~
>
> The output should be a job ID.
> With the job ID you can query the status of your "job" on the cluster
>
> ~~~
> $ qstat 125932.pbs02 
> ~~~
> 
> If all went well, the output of `qstat` should show `R` (running)
> under that `S` (status) column:
>
> ~~~
> Job id            Name             User              Time Use S Queue
> ----------------  ---------------- ----------------  -------- - -----
> 1453338.pbs02     sleep            atrikut           00:00:00 R c1_solo>
> ~~~
>
> Eventually, your job will complete, yielding the following messge:
>
> ~~~
> qstat: 1453338.pbs02 Job has finished, use -x or -H to obtain historical job information
> ~~~
>
> Once a job is completed, you can check for additional details using:
>
> ~~~
> qstat -xf <jobID>
> ~~~
>
> What node did your jobs run on? How much walltime did they take?
> Does this align with the expected walltime for this job?
{: .challenge}

Once a batch job completes,
you should also see two new files in the directory from which
`qsub` was run:

~~~
$ ls 

sleep.e1453338  sleep.o1453338  ....
~~~

These two files (`<job_name>.o<job_ID> and <job_name>.e<job_ID>`) contain
the standard output (i.e., any text printed as output from the commands),
and the standard error (i.e., any error text printed from the commands)
of the job respectively.

First, let's check that no errors were produced:

~~~
cat sleep.e1453338
~~~

And let's also have a look at the output:

~~~
Good morning!!


+------------------------------------------+
| PALMETTO CLUSTER PBS RESOURCES REQUESTED |
+------------------------------------------+

mem=1gb,ncpus=1,walltime=00:05:00


+-------------------------------------+
| PALMETTO CLUSTER PBS RESOURCES USED |
+-------------------------------------+

cpupercent=0,cput=00:00:00,mem=5456kb,ncpus=1,vmem=338432kb,walltime=00:01:41

~~~

We see that a job summary is printed in addition to the standard output.
