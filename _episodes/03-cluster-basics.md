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
