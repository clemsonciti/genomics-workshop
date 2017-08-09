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



