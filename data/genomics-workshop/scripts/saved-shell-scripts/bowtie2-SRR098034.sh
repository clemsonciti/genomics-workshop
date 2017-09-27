#!/bin/bash

bowtie2 -x Reference/Ecoli_Ref.fa \
	--rg-id SRR098034 --rg "SM:ZDB83" \
	-U Trimming/SRR098034.trim.fq \
	-S Alignment/SRR098034.sam
