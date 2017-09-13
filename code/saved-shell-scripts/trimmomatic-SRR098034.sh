#!/bin/bash

java -jar /software/trimmomatic/0.30/trimmomatic-0.30.jar SE \
	-threads 2 \
	-phred33 \
	-trimlog Trimming/SRR098034.log \
	Raw_Fastq/SRR098034.fastq Trimming/SRR098034.trim.fq \
	ILLUMINACLIP:adapters.fasta:2:30:10 \
	TRAILING:20 MINLEN:30 AVGQUAL:30
