---
title: "Parallel Variant Calling Workflow"
teaching: 20
exercises: 1
---

# Overview
The purpose of today's exercises will be to perform all of the steps in the Variant Calling Workflow on the 4 remaining raw data files.  At each step, we will look at a different possible way for parallelizing a process to run on all 4 files simultaneously.  After running the GATK HaplotypeCaller on all files, we will run the final step of the variant calling pipeline (GATK GenotypeGVCFs) in order to generate the final file containing polymorphism information for each individual.  As a reminder, the steps in the variant calling workflow are:

1. Quality check with FastQC
2. Trim and filter with Trimmomatic
3. Mapping with Bowtie2
4. Format and Sort with Samtools 
5. Variant calling with GATK:
	1. Get genotypes in GVCF format for each individual (GATK HaplotypeCaller)
	2. Perform joint genotyping across all .gvcf files to get the final VCF file for the population (GATK GenotypeGVCFs)

# Before Starting
As part of the exercises from Day 1, you should now have all 5 fastq files in your Raw_Fastq folder:
~~~
cd ~/genomics-workshop/Raw_Fastq
ls -lh
[ecoope4@node0071 Raw_Fastq]$ ls -lh
total 2.5G
-rwxr-xr-x 1 ecoope4 cuuser 1.5G Aug  8  2017 SRR098034.fastq
-rwxr-xr-x 1 ecoope4 cuuser 1.5G Jul 18 14:43 SRR098035.fastq
-rwxr-xr-x 1 ecoope4 cuuser 1.4G Jul 18 14:43 SRR098038.fastq
-rwxr-xr-x 1 ecoope4 cuuser 1.6G Jul 18 14:44 SRR098039.fastq
-rwxr-xr-x 1 ecoope4 cuuser 1.5G Jul 18 14:45 SRR098289.fastq 
~~~

If you do NOT have all of the fastq file successfully downloaded, then you can get a copy of them here:
(Put link to zipped folder of sample files!!!)

To unzip the folder and put the files in the correct directory (only necessary if you don't already have the files), do:
~~~
unzip SampleFiles.zip
cd SampleFiles/
mv * ~/genomics-workshop/Raw_Fastq/
~~~

# Trimmomatic (Using a For Loop)
For the purposes of this workshop, we are going to skip the first step of running FastQC on each file, and instead proceed directly to Step 2: Trimmomatic.

In this example, we're going to use a For loop to run the same Trimmomatic command on every file inside of the Raw_Fastq directory.  The loop will help us do this all with one script (without having to write the command out 4 separate times), but note that in this first instance we are NOT running anything in parallel!  The command will run sequentially on each file (one after the other), but only one command will be executed at a time.  For the case of Trimmomatic, this is still fairly efficient, because Trimmomatic itself was written to be able to use multiple cores in one run (with the -threads option), but be aware that this will not be the case for all programs.

First, we will submit the PBS job script for running Trimmomatic (it will take about 15-17 minutes to run completely on all 5 files), and while it is running, we will look at the contents of the script more closely.

(NOTE FOR LIZ AND ASHWIN: WE NEED TO GO IN AND CHANGE THE PATHS IN THIS SCRIPT BEFORE INCLUDING IN THE FINAL WORKSHOP FOLDER)

~~~
[ecoope4@login001 scripts]$ qsub trimmomatic.sh 
890699.pbs02
~~~

Once you have the script running, use the "more" command to open up the script and look at it:
~~~
130-127-150-127:code lizcooper$ more trimmomatic.sh 
#!/bin/bash
#PBS -N Trimm
#PBS -l select=1:ncpus=8:mem=11gb,walltime=2:00:00
#PBS -j oe


echo "START ------------------------------"

module add java/1.8.0

src=/zfs/tillers/liz/workshop
export adapt=$src/adapters.fasta

for file in $src/Raw_Fastq/*.fastq
do
        export prefix=`basename $file .fastq`
        export outname="$prefix.trim.fq"
        java -jar /software/trimmomatic/0.30/trimmomatic-0.30.jar SE \
                -threads 8 \
                -phred33 \
                -trimlog "$src/Trimming/$prefix.log" \
                $file $src/Trimming/$outname \
                ILLUMINACLIP:$adapt:2:30:10 \
                TRAILING:20 MINLEN:30 AVGQUAL:30 
done    

echo "FINISH ----------------------------"
~~~

