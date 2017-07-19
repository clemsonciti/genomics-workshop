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
[nelle@login001 genomics-workshop]$ mkdir Raw_Fastq
~~~

Log onto an interactive node:
~~~
[nelle@login001 genomics-workshop]$ qsub -I -l select=1:ncpus=2:mem=31gb:interconnect=1g,walltime=1:00:00
qsub (Warning): Interactive jobs will be treated as not rerunnable
qsub: waiting for job 465026.pbs02 to start
qsub: job 465026.pbs02 ready

[nelle@node1134 ~]$
~~~

Load the SRA toolkit module:
~~~ 
[nelle@node1134 ~]$ module load sratoolkit/2.8.2-1
~~~

Now, fetch the first .sra file from the ncbi repository:
~~~
[nelle@node1134 ~]$ prefetch -v SRR098034

2017-07-18T13:01:54 prefetch.2.8.2: KClientHttpOpen - connected to www.ncbi.nlm.nih.gov
2017-07-18T13:01:54 prefetch.2.8.2: KClientHttpOpen - verifying CA cert
2017-07-18T13:01:54 prefetch.2.8.2: KClientHttpOpen - connected to sra-download.ncbi.nlm.nih.gov
2017-07-18T13:01:54 prefetch.2.8.2: KClientHttpOpen - verifying CA cert
2017-07-18T13:01:54 prefetch.2.8.2: 1) Downloading 'SRR098034'...
2017-07-18T13:01:54 prefetch.2.8.2:  Downloading via https...
2017-07-18T13:01:54 prefetch.2.8.2: https://sra-download.ncbi.nlm.nih.gov/traces/sra1/SRR/000095/SRR098034 -> /home/nelle/ncbi/public/sra/SRR098034.sra.tmp.22825.tmp
2017-07-18T13:01:54 prefetch.2.8.2: KClientHttpOpen - connected to sra-download.ncbi.nlm.nih.gov
2017-07-18T13:01:54 prefetch.2.8.2: KClientHttpOpen - verifying CA cert
2017-07-18T13:01:59 prefetch.2.8.2: /home/nelle/ncbi/public/sra/SRR098034.sra.tmp.22825.tmp (0)
2017-07-18T13:01:59 prefetch.2.8.2: 1) 'SRR098034' was downloaded successfully
[nelle@node1134 ~]$
~~~

The .sra file is automatically put in a specific place created in your home directory, which can be found at:
~~~
[nelle@node1134 ~]$ cd ~/ncbi/public/sra/
[nelle@node1134 ~]$ ls
SRR098034.sra
~~~

From this directory, convert the file to .fastq format, and have it output to the Raw data folder we made inside of our workshop directory:
~~~
[nelle@node1134 sra]$ fastq-dump SRR098034.sra --outdir ~/genomics-workshop/Raw_Fastq/
Read 7575758 spots for SRR098034.sra
Written 7575758 spots for SRR098034.sra
~~~

This step might take a few minutes to complete.  Once it is finished, cd into the Raw_Fastq directory and check that the output file is there.

~~~
[nelle@node1134 sra]$ cd ~/genomics-workshop/Raw_Fastq/
[nelle@node1134 Raw_Fastq]$ ls
SRR098034.fastq
~~~

~~~
[nelle@node1134 Raw_Fastq]$ head SRR098034.fastq
@SRR098034.1 30KR6AAXXLesnki_set_1_2:1:2:1394:602 length=36
GGTTGACGAGTTTATCGAAGTTTTTCGCCAGAACCA
+SRR098034.1 30KR6AAXXLesnki_set_1_2:1:2:1394:602 length=36
CCCCCCCCCCCCCCCCCCCCCCCCCC9CAC>??.A>
@SRR098034.2 30KR6AAXXLesnki_set_1_2:1:2:1283:96 length=36
GATAAGCCCGATTTGAAGGCATAGTTTACCATGCGC
+SRR098034.2 30KR6AAXXLesnki_set_1_2:1:2:1283:96 length=36
CCCCCCCA6ACCCCACCAC?CCC<CCCCA<AC?<?%
@SRR098034.3 30KR6AAXXLesnki_set_1_2:1:2:1255:120 length=36
GTTCGCCAAGTTCATTGGCAGAGGGAGCCAACAGAC
~~~
