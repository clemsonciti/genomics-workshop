---
title: "Variant calling workflow"
teaching: 20
exercises: 1
---

# Overview
A variant calling workflow starts with raw sequencing data for multiple individuals, and ends with a single file containing only the genomic positions where at least one individual in the population has a mutation or polymorphism.  The steps in this process are:

1. Using FastQC to verify the quality of each raw data file.
2. Removing low quality and spurious sequences (Trimmomatic).
3. Mapping individual sequences to a reference genome (Bowtie2)
4. Convert SAM files produces in Step 3 to BAM format files, then sort and index them (Samtools)
5. Perform variant calling on the sorted BAM files (GATK)
	1. Get genotypes in GVCF format for each individual (GATK HaplotypeCaller)
	2. Perform joint genotyping across all .gvcf files to get the final VCF file for the population (GATK GenotypeGVCFs)

On Day 1, we will only be working with a single file to start with, so we will *not* do Step 5, part 2 (yet).  On Day 2, after we have looked at running those same steps in parallel on all of the other raw data files, we will complete the last part of the pipeline to get the final VCF file.

# FastQC
To use the program FastQC, we need to make sure that we are logged onto palmetto with X11 forwarding (or the Windows equivalent):
~~~
ssh -X ecoope4@user.palmetto.clemson.edu
~~~

We will also need to use the -X option when requesting an interactive node:
~~~
qsub -X -I -l select=1:ncpus=2:mem=31gb:interconnect=1g,walltime=1:00:00
~~~

Once you are on an interactive node, load the software module for FastQC, and start the program by simply typing "fastqc" in the command line:
~~~
module load fastqc/0.10.1
fastqc
~~~

This will open up an X11 window with the FastQC startup window.  From the File menu, select Open, and then navigate to the SRR098034.fastq file in your genomics-workshop/Raw_Fastq folder.

FastQC checks multiple aspects of the raw data: (list these later)

# Trimming
Trimmomatic will clean up the raw sequencing data by removing any reads or parts of reads that fall below a specified quality threshold, and also by checking for and removing adapter sequences left over from the Illumina library construction.

Trimmomatic is a Java application, so we need to first load the java module.  When we run the program, note that we need to specify the path to the .jar file:

~~~
cd ~/genomics-workshop
mkdir Trimming

module load java/1.8.0

java -jar /software/trimmomatic/0.30/trimmomatic-0.30.jar SE \
-threads 2 \
-phred33 \
-trimlog Trimming/SRR098034.log \
Raw_Fastq/SRR098034.fastq Trimming/SRR098034.trim.fq \
ILLUMINACLIP:adapters.fasta:2:30:10 \
TRAILING:20 MINLEN:30 AVGQUAL:30
~~~

Need to go back and add some notes on the parameters...

# Mapping to the Reference Genome
For the sample data used in this workshop, we will be mapping everything to the genome of *E.coli* strain REL606.  The fasta file for this genome has already been provided to you, but it was originally downloaded from [Ensembl](http://www.ensembl.org/index.html), which is an excellent resource for many publicly available genomes.  

We will be using the program Bowtie2 for mapping, which is a very commonly used aligner.  It is also very similar to BWA.  Bowtie2 will produce SAM format files, which has become the standard output format for mapped reads, and is now used by many programs.

Before we can align any reads, we have to first create an Index for our reference genome.  You only need to do this once for any reference genome file (i.e. not every time you map to it).
~~~
cd ~/genomics-workshop
module load bowtie2/2.1.0
bowtie2-build Reference/Ecoli_Ref.fa Reference/Ecoli_Ref
~~~

Now that we have prepared the reference genome, we can run the alignment step on our trimmed .fastq file (this will take ~8 minutes):
~~~
mkdir Alignment

bowtie2 -x Reference/Ecoli_Ref \
--rg-id SRR098034 --rg "SM:ZDB83" \
-U Trimming/SRR098034.trim.fq \
-S Alignment/SRR098034.sam
~~~

# Samtools
After creating a SAM file with an alignment program, there can often be one or more intermediate steps you need to run before you can use the file in another program.  In practice, the exact steps you need/want will depend on your data set and what you want to do with it, but the [Samtools](http://samtools.sourceforge.net/) toolkit contains numerous useful functions specifically for SAM files.  The functions we will use today are:
1. Converting from SAM to BAM format (BAM is just a compressed form of SAM)
2. Sorting the BAM file by numerical position
3. Creating an index for the BAM file, which is required by GATK
4. Creating 2 new index files for the reference (in formats used by Samtools and GATK)

The command lines for running these 4 steps are:
~~~
mkdir BamFiles
module load samtools/1.4

samtools view -O BAM Alignment/SRR098034.sam -o BamFiles/SRR098034.bam
samtools sort -o BamFiles/SRR098034.sort.bam -O BAM BamFiles/SRR098034.bam
samtools index -b BamFiles/SRR098034.sort.bam
samtools faidx Reference/Ecoli_Ref.fa
samtools dict -o Reference/Ecoli_Ref.dict Reference/Ecoli_Ref.fa
~~~

# GATK HaplotypeCaller
The final task for Day 1 is to run the GATK HaplotypeCaller to generate a .gvcf file for our first individual. This file will be used on Day 2, once we have also created the .gvcf files for the 4 other individuals, to generate the final VCF file.

[GATK](https://software.broadinstitute.org/gatk/), like Samtools, has a lot of uses and functions in additions to the ones we will use in this workshop.  It is one of the most commonly used programs for variant calling (even though there are many others), which is why we are using it here.  To run the HaplotypeCaller:

~~~
module load java/1.8.0
mkdir SNPs

java -jar GenomeAnalysisTK.jar -T HaplotypeCaller \
-I BamFiles/SRR098034.sort.bam \
-R Reference/Ecoli_Ref.fa \
-o SNPs/SRR098034.gvcf \
-ERC GVCF \
-variant_index_type LINEAR \
-variant_index_parameter 128000 \
-ploidy 1
~~~

Need to add additional notes on parameters...

