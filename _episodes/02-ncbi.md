---
title: "Navigating the NCBI database"
teaching: 20
exercises: 1
questions:
- "How do I navigate the NCBI database?"
- "How do I download an SRA file?"
keypoints:
- "Keypoint 1"
- "Keypoint 2"
---

# Data Set
In this workshop, we will be working with *E. coli* data generated as part of a [2012 study](https://www.nature.com/nature/journal/v489/n7417/full/nature11514.html) that tracked the evolutionary origins of a novel trait.  For our purposes, we'll only be using 5 of the sequencing libraries, described in the table below.

| SRA ID | Sample Name | Sequencing Description |
|--------|:-----------:|:----------------------:|
| SRR098034 | ZDB83 | Illumina Single End |
| SRR098035 | ZDB87 | Illumina Single End |
| SRR098038 | ZDB107 | Illumina Single End |
| SRR098039 | ZDB111 | Illumina Single End |
| SRR098289 | ZDB564 | Illumina Single End |

# Finding the Data Online

# Download the Data
First, create a directory to hold the raw data files:
~~~
mkdir Raw_Fastq
~~~

Log onto an interactive node:
~~~
qsub -I -l select=1:ncpus=2:mem=31gb:interconnect=1g,walltime=1:00:00
~~~

Load the SRA toolkit module:
~~~ 
module load module load sratoolkit/2.8.2-1
~~~

Now, fetch the first .sra file from the ncbi repository:
~~~
prefetch -v SRR098034 
~~~

The .sra file is automatically put in a specific place created in your home directory, which can be found at:
~~~
cd ~/ncbi/public/sra
ls
~~~

From this directory, convert the file to .fastq format, and have it output to the Raw data folder we made inside of our workshop directory:
~~~
fastq-dump SRR098034.sra --outdir ~/genomics-workshop/Raw_Fastq/
~~~

This step might take a few minutes to complete.  Once it is finished, cd into the Raw_Fastq directory and check that the output file is there.
  
