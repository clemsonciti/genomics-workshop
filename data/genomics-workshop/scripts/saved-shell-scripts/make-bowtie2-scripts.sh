#!/bin/bash

export x=1
while read -r line
do
	export sraID=`echo $line | cut -d" " -f1`
	export sample=`echo $line | cut -d" " -f2`
	sed "s/TEMP_SRA/$sraID/g" ../scripts/bowtie2-aln.sh >scripts/temp.sh
	sed "s/TEMP_SAMPLE/$sample/g" ../scripts/temp.sh >"scripts/bowtie2.$x.sh"
	let x=(x+1)
done < data/fileList.txt 
