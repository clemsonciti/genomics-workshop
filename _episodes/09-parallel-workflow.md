---
title: "Parallel Variant Calling Workflow"
teaching: 20
exercises: 1
---

# Overview
The purpose of these exercises will be to perform all of the steps in the Variant Calling Workflow on the 4 remaining raw data files.  At each step, we will look at a different possible way for parallelizing a process to run on all 4 files simultaneously.  After running the GATK HaplotypeCaller on all files, we will run the final step of the variant calling pipeline (GATK GenotypeGVCFs) in order to generate the final file containing polymorphism information for each individual.  As a reminder, the steps in the variant calling workflow are:

1. Quality check with FastQC
2. Trim and filter with Trimmomatic
3. Mapping with Bowtie2
4. Format and Sort with Samtools 
5. Variant calling with GATK:
	1. Get genotypes in GVCF format for each individual (GATK HaplotypeCaller)
	2. Perform joint genotyping across all .gvcf files to get the final VCF file for the population (GATK GenotypeGVCFs)

# Before Starting
As a reminder, you should now have all 5 fastq files in your Raw_Fastq folder:
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

In this example, we'd like to use a `for` loop to run the same Trimmomatic command on every file inside of the Raw_Fastq directory.  The loop will help us do this all with one script (without having to write the command out 4 separate times), but note that in this first instance we are NOT running anything in parallel!  The command will run sequentially on each file (one after the other), but only one command will be executed at a time.  For the case of Trimmomatic, this is still fairly efficient, because Trimmomatic itself was written to be able to use multiple cores in one run (with the `-threads` option), but be aware that this will not be the case for all programs.

 Since you have some experience with for loops from yesterday, try drafting a PBS script with a for loop that runs the trimmomatic command we used earlier on all the files. Don't for get to add the necessary #PBS headings at the top and we recommend selecting 8 cpus, 11gb of memory, and a walltime of 2 hours. Once you have a script drafted try to submit it and see if it runs.

Now you can compare to the script below which you can find in the scripts directory
~~~
130-127-150-127:code lizcooper$ more trimmomatic.sh 
#!/bin/bash
#PBS -N Trimm
#PBS -l select=1:ncpus=8:mem=11gb,walltime=2:00:00
#PBS -j oe


echo "START ------------------------------"

module add java/1.8.0

src=/home/$USER/genomics-workshop
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

The Trimmomatic command itself is relatively unchanged from the command we ran in interactive mode (all of the parameter options are the same), but there are a few key pieces of this script that help with running the command effectively in a loop.  First, note these two lines:
~~~
src=/home/$USER/genomics-workshop
export adapt=$src/adapters.fasta
~~~

The first line saves a filepath to the variable "$src."  The second line uses that "$src" variable, and also points to the file with the adapters in it to save that whole file path as the variable "$adapt."  Both of these lines are optional (we could just call up the full file path every time we refer to any file), but you should notice how they help save space by avoiding having to type long file names later in the script.

The next key piece:
~~~
for file in $src/Raw_Fastq/*.fastq
~~~
tells the loop to operate on EVERY file inside of the Raw_Fastq folder that ends with the ".fastq" extension.  Note the use of the "*" wildcard character, which will match anything coming before ".fastq."

In order to automatically get a name for each of my output files that will correspond to the names of each input file, I first use the `basename` command, which will take off all of the extra path information before the file name itself.  Here, I am also specifying a suffix (.fastq), so that will be taken off by `basename` as well.  This will leave me with just my sample name, which I can use as a prefix for my output file name in the next line:
~~~
export prefix=`basename $file .fastq`
export outname="$prefix.trim.fq"
~~~

Finally, Trimmomatic is run with the same command line we used before, but now we have the variables "$file" specifying the input, "$outname" specifying the output, and "$adapt" specifying the adapters.fasta file. 

# Run Bowtie2 with multiple job scripts
Instead of writing a single loop that will process all files sequentially, you can take advantage of the many nodes available on Palmetto and submit a separate job for each file.  These will run simultaneously, thereby finishing faster than if they ran one after the other.  Note that you want to make sure you are requesting the appropriate resources for each individual job (and not taking up nodes that you don't need!).

One potential drawback of this approach is that you actually have to write a separate PBS script for each input file.  This doesn't seem too bad for 4 files, but could become very time consuming for a large number of files.  To show you how you can automate this process, we will start with a template PBS script for running Bowtie2, and then we will run **another** shell script (just a loop) to search and replace filenames in order to make a new script for each input file.  The template script is `bowtie2-aln.sh`, so open it and look at how it is written:
~~~
[ecoope4@node0050 scripts]$ less bowtie2-aln.sh 
#!/bin/bash
#PBS -N bowAln
#PBS -l select=1:ncpus=8:mem=11gb,walltime=2:00:00
#PBS -j oe


echo "START ------------------------------"

module add bowtie2/2.1.0 

src=/home/$USER/genomics-workshop
export srrname=TEMP_SRA
export sample=TEMP_SAMPLE

bowtie2 -p 8 -x $src/Reference/E_coli \
        --rg-id  $srrname \
        --rg "SM:$sample" \
        -U "$src/Trimming/$srrname.trim.fq" \
        -S "$src/Alignment/$srrname.sam"


echo "FINISH ----------------------------"
~~~

In this script, you can see that we define 2 variables, $srrname and $sample, that are each given a "TEMP" value to start with.  These are the values we will replace to be different for each input file.  Since the rest of the script uses the variable names, these will point to the correct values every time.  The parameters for Bowtie2 are the same as we used earlier

The 2 values that we want to replace in this script are the SRA database ID for each file, and the Sample Name that corresponds to each ID.  These are given in a table describing the data, and are also provided in a text file called fileList.txt:
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
do
	export sraID=`echo $line | cut -d" " -f1`
	export sample=`echo $line | cut -d" " -f2`
	sed "s/TEMP_SRA/$sraID/g" ../scripts/bowtie2-aln.sh >../scripts/temp.sh
	sed "s/TEMP_SAMPLE/$sample/g" ../scripts/temp.sh >"../scripts/bowtie2.$x.sh"
	let x=(x+1)
	done <fileList.txt 
~~~

When this runs, you should have 5 new bowtie2 scripts.  Try using less to look at a couple of them, and check that the SRA ids and sample names have been changed correctly in each of them.  Now, you can submit each of these (you actually do NOT need to submit `bowtie2.1.sh`, since we ran this interactively) using the `qsub` command.

Despite how cumbersome this method may seem at first, this is an effective way to deal with a lot of files.  This method works the same no matter how many individual files you have, and having a separate job run for each file makes it easy to keep track of errors if anything goes wrong with just one or a few files. 

#  Using a single loop and job control to run Samtools

In both the Trimmomatic and Bowtie2 scripts we used above, we were able to take advantage of the fact that those two programs have multithreading options (`-threads` in Trimmomatic and `-p` in Bowtie2) in order to use all of the available cores on a given node.  However, sometimes you will need to use programs that don't permit multithreading, so it is worth considering alternate strategies for running these efficiently on the cluster.  

The unix idea of "job control" has options that can allow a user to send a running process to the "background," and then start a new task before the first one has finished.  This means that you can run multiple processes simultaneously, even on a single node.

In this example, we are going to use a loop that will START a Samtools command on each SAM file in our Alignment directory, and then send each process to the background to run while the next file is started.  In this way, Samtools will be running on all of our files simultaneously, and the script will finish once all of the background processes have been completed.  

(As a side note: the current version of Samtools actually does have a multithreading option, even though we are not using it here.  Previous versions of Samtools, however, did NOT have multithreading capability, so the use of job control used to be crucial with Samtools.)

The script for running Samtools as background processes is `samtools.sh.`  Use the `more` command to look at the contents of this script.

~~~
$ more samtools.sh 
#!/bin/bash
#PBS -N sam
#PBS -l select=1:ncpus=8:mem=11gb,walltime=2:00:00
#PBS -j oe


echo "START ------------------------------"

module add samtools/1.4 

src=/home/$USER/genomics-workshop

### Loop over all aligned SAM files
### Convert to bam and sort
for file in $src/Alignment/*.sam
do
	### Get the sample prefix from the file name
	export sample=`basename $file .sam`
	export outname="$sample.sort.bam"
	samtools view -O BAM $file | samtools sort -o $src/BamFiles/$outname -O BAM - &
done

wait

### This is a quick loop to get the index files needed for GATK
for file in $src/BamFiles/*.sort.bam
do
	samtools index -b $file &
done
wait

echo "FINISH ----------------------------"
~~~

The line doing all of the work in this script is:
~~~
samtools view -O BAM $file | samtools sort -o $src/BamFiles/$outname -O BAM - &
~~~

The `&` at the end of the line is what tells the program to send that process into the background (and thus start on the next process without waiting on this one to finish).  You should also note the `wait` command that comes right after the loop; this tells the program to make sure all of the background processes are finished before going on to the next loop (and before exiting the script).

This samtools command line is also different from the one we saw on Day 1 in that it is using a `|` to combine the "view" and "sort" commands into 1 line.  This isn't necessary for running the command in the background; it is just an extra trick to show you how pipes can be useful to simplify your code and avoid creating uneeded intermediate files.

Let's look closer at how the `|` command line is working.  First, recall that the samtools "view" command line we used on Day 1 looked like this:
~~~
samtools view -O BAM <InFile> -o <OutFile>
~~~
In the first half of the piped command, we no longer have the `-o <OutFile>` option for the "view" command, which means that by default the output of samtools would go to STDOUT, instead of into a file.  The `|` is taking this output from STDOUT and redirecting it straight into STDIN for the "sort" command.

The "sort" command that we used on Day 1 looked like this:
~~~
samtools sort -o <OutFile> -O BAM <InFile>
~~~

In our piped command line, we are still writing the final output to a file, so the `-o` option is still used.  The change is that we want the input to come from STDIN instead of from a file, BUT we also have to deal with the fact the Samtools requires an input file name as an argument, so we can't just leave this blank.  The `-` character at the end of the line is acting as a special placeholder for the filename argument that directs the program to use STDIN instead.

After looking over how `samtools.sh` is written, submit the job using the `qsub` command to get sorted BAM files and BAM index (.bai) files for every sample.  This should take about 4-5 mintues to run.

#  Using Gnu-Parallel to call haplotypes with GATK
The last type of parallelization technique that we will look at is using GNU Parallel.  [GNU-parallel](https://www.gnu.org/software/parallel/) is a separate software written for the purpose of simplifying the parallel execution of command line jobs and scripts. Gnu-parallel may seem intimidating, but is one of the most effective ways of parallelizing tasks on the palmetto cluster. There are several different ways that one can call gnu-parallel with slightly different command structures. In my opinion, one of the most straightforward methods is to piping to direct gnu-parallel. To use this method it is best to start out by creating a list file that passes each item in the list to gnu parallel. For instance, if we wanted to run a command in parallel on three files called red, blue, and yellow, we could navigate to th directory containing them and do this:
~~~
[user@node]$cd colors
[user@node]$ls > list.txt
~~~
If we were to cat that list we should see:
~~~
blue
red
yellow
~~~
Now, if we wanted to execute a task like echoing the color we could pipe that output into parallel in the following way:

~~~
[usr@node]cat list.txt | parallel "echo {}"
~~~

With this format of parallel, the commands inside the quotes are what is passed to the machine and the bracket is filled with the output from the pipe. You can ask parallel to show you this if you specify it using the -v or --verbose switch, which should look like this:

~~~
[usr@node]cat list.txt | parallel --verbose "echo {}"
echo blue
echo red
echo yellow
blue
red
yellow
~~~

Using the verbose option of parallel is a great way to trouble shoot, so I recommend using it most of the time. Additionally, because the brackets are replaced exactly with what is piped in it is very important to make the sure the list contents are exactly what you want. For instance, if you would like to use a file blue.txt as input and have it output as blue.out it would be easier to have the list contents be "blue" instead of "blue.txt". There are a variety of strategies to accomplish this, especially using piping and editing commands like sed. 

The way that we will be implementing gnu-parallel in this exercise is to run GATK in parallel over a list of files. Today's example is intended to give you a brief introduction of its syntax and how to use it; the above link to the homepage provides access to more detailed tutorials if you want to know more about different ways to implement gnu-parallel.  

The script we will use for running the GATK HaplotypeCaller on all sorted BAM files is `gatk-hapCall.sh`, so let's first look at its contents:
~~~
less gatk-hapCall.sh

#!/bin/bash
#PBS -N gatk
#PBS -l select=1:ncpus=8:mem=15gb,walltime=2:00:00
#PBS -j oe


echo "START ------------------------------"

module add java/1.8.0 
module load gnu-parallel
module load GATK

src=/home/$USER/genomics-workshop

### Use gnu-parallel to use multiple cores
### within one script
cd $src/BamFiles
cat list.txt | parallel --verbose -j 8 "java -jar $GATK \
        -T HaplotypeCaller \
        -R $src/Reference/Ecoli_Ref.fa \
	-ploidy 1 \
        -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 \
        -I {}.sort.bam -o $src/SNPs/{}.gvcf" 

echo "FINISH ----------------------------"
~~~

There are a few key pieces of this script to take note of.  First, the line:
~~~
module load gnu-parallel
~~~
is loading the software module that has been pre-installed on Palmetto.  The java module is also loaded because it is required by GATK.  The next thing to notice is how the program is called: `parallel --plus -j 8` is the command that calls the gnu-parallel program.  The `-j 8` option tells the program we want to run 8 jobs at a time in parallel.  

The parallel options get called just before we call the actual GATK tool (with the `java -jar` syntax that we used on Day 1).  The rest of the GATK command line stays mostly the same. The last thing that we need is the list generated in the correct directory.

Based on the structure of the parallel command try to generate the list.txt file that we will need in th directory we need it in. Try to avoid typing out the sample names if you can and instead use some of the commandline tricks we have learned so far. Once you think you have got it submit the `gatk-hapCall.sh` job script to try and get the .gvcf files for each individual. If it doesn't work immediately try to trouble shoot.

# GATK GenotypeGVCFs
The very last step in the pipeline, now that we have multiple .gvcf files, is to generate a single SNP file containing the genotypes at all sites with a sequence variant (mutation) in one or more individuals in our sample.  This is done with the `GATK` tool `GenotypeGVCFs`:

~~~

java -jar $GATK \
	-T GenotypeGVCFs \
	-R Reference/Ecoli_Ref.fa \
	--variant SNPs/SRR098034.gvcf \
	--variant SNPs/SRR098035.gvcf \
	--variant SNPs/SRR098038.gvcf \
	--variant SNPs/SRR098039.gvcf \
	--variant SNPs/SRR098289.gvcf \
	-o SNPs/combined.vcf
~~~

Note that in the above command line each individual file you want to include is specified with the `--variant` flag, and you can include as many individuals as you want in a single VCF file.

Create your own PBS job script to run this final step (make sure you specify your paths correctly), and then submit it to get the final SNP file.  This should take 10 minutes (or less) to run.
