#!/bin/bash
#PBS -N bowAln
#PBS -l select=1:ncpus=8:mem=11gb,walltime=2:00:00
#PBS -j oe


echo "START ------------------------------"

module add bowtie2/2.1.0 

src=/scratch2/$USER/genomics-workshop
export srrname=TEMP_SRA
export sample=TEMP_SAMPLE

bowtie2 -p 8 -x $src/Reference/Ecoli_Ref \
	--rg-id  $srrname \
	--rg "SM:$sample" \
	-U "$src/Trimming/$srrname.trim.fq" \
	-S "$src/Alignment/$srrname.sam"

echo "FINISH ----------------------------"



