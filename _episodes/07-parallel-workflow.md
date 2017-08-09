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

The Trimmomatic command itself is relatively unchanged from the command we ran on Day 1 in interactive mode (all of the parameter options are the same), but there are a few key pieces of this script that help with running the command effectively in a loop.  First, note these two lines:
~~~
src=/zfs/tillers/liz/workshop
export adapt=$src/adapters.fasta
~~~

The first line saves a filepath to the variable "$src."  The second line uses that "$src" variable, and also points to the file with the adapters in it to save that whole file path as the variable "$adapt."  Both of these lines are optional (we could just call up the full file path every time we refer to any file), but you should notice how they help save space by avoiding having to type long file names later in the script.

The next key piece:
~~~
for file in $src/Raw_Fastq/*.fastq
~~~
tells the loop to operate on EVERY file inside of the Raw_Fastq folder that ends with the ".fastq" extension.  Note the use of the "*" wildcard character, which will match anything coming before ".fastq."

In order to automatically get a name for each of my output files that will correspond to the names of each input file, I first use the "basename" command, which will take off all of the extra path information before the file name itself.  Here, I am also specifying a suffix (.fastq), so that will be taken off by basename as well.  This will leave me with just my sample name, which I can use as a prefix for my output file name in the next line:
~~~
export prefix=`basename $file .fastq`
export outname="$prefix.trim.fq"
~~~

Finally, Trimmomatic is run with the same command line we used before, but now we have the variables "$file" specifying the input, "$outname" specifying the output, and "$adapt" specifying the adapters.fasta file. 

# Bowtie2 (with multiple job scripts)
Instead of writing a single loop that will process all files sequentially, you can take advantage of the many nodes available on Palmetto and submit a separate job for each file.  These will run simultaneously, thereby finishing faster than if they ran one after the other.  Note that you want to make sure you are requesting the appropriate resources for each individual job (and not taking up nodes that you don't need!).

One potential drawback of this approach is that you actually have to write a separate PBS script for each input file.  This doesn't seem too bad for 4 files, but could become very time consuming for a large number of files.  To show you how you can automate this process, we will start with a template PBS script for running Bowtie2, and then we will run another shell script to search and replace filenames in order to make a new script for each input file.  The template script is bowtie2-aln.sh, so open it and look at how it is written:
~~~
[ecoope4@node0050 scripts]$ less bowtie2-aln.sh 
#!/bin/bash
#PBS -N bowAln
#PBS -l select=1:ncpus=8:mem=11gb,walltime=2:00:00
#PBS -j oe


echo "START ------------------------------"

module add bowtie2/2.1.0 

src=/zfs/tillers/liz/workshop
export srrname=TEMP_SRA
export sample=TEMP_SAMPLE

bowtie2 -p 8 -x $src/Reference/E_coli \
        --rg-id  $srrname \
        --rg "SM:$sample" \
        -U "$src/Trimming/$srrname.trim.fq" \
        -S "$src/Alignment/$srrname.sam"


echo "FINISH ----------------------------"
~~~

In this script, you can see that we define 2 variables, $srrname and $sample, that are each given a "TEMP" value to start with.  These are the values we will replace to be different for each input file.  Since the rest of the script uses the variable names, these will point to the correct values every time.  The parameters for Bowtie2 are the same as we used on Day 1.

The 2 values that we want to replace in this script are the SRA database ID for each file, and the Sample Name that corresponds to each ID.  These are given in a table from Day 1 describing the data, and are also provided in a text file called fileList.txt:
~~~
[ecoope4@login001 workshop]$ less fileList.txt 
SRR098034       ZDB83
SRR098035       ZDB87
SRR098038       ZDB107
SRR098039       ZDB111
SRR098289       ZDB564
~~~

We are going to read in each line of this file, and then make a new PBS script with the corresponding SRA and Sample IDs for each line.  I also want to make sure that each new script that gets written has a new name (so I don't overwrite the same script over and over again).  We'll use a numeric "counter" variable, $x, to help us do this.  The loop for writing the new scripts is:
~~~
[ecoope4@node0050 workshop]$ export x=1
[ecoope4@node0050 workshop]$ while read -r line
> do
> export sraID=`echo $line | cut -d" " -f1`
> export sample=`echo $line | cut -d" " -f2`
> sed "s/TEMP_SRA/$sraID/g" scripts/bowtie2-aln.sh >scripts/temp.sh
> sed "s/TEMP_SAMPLE/$sample/g" scripts/temp.sh >"scripts/bowtie2.$x.sh"
> let x=(x+1)
> done <fileList.txt 
~~~

When this runs, you should have 5 new bowtie2 scripts.  Try using less to look at a couple of them, and check that the SRA ids and sample names have been changed correctly in each of them.  Now, you can submit each of these (you actually do NOT need to submit bowtie2.1.sh, since we ran this file yesterday) using the qsub command. 
