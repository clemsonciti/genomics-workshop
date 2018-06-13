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



