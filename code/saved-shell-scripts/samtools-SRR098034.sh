#!/bin/bash

samtools view -O BAM Alignment/SRR098034.sam -o BamFiles/SRR098034.bam
samtools sort -o BamFiles/SRR098034.sort.bam -O BAM BamFiles/SRR098034.bam
samtools index -b BamFiles/SRR098034.sort.bam
samtools faidx Reference/Ecoli_Ref.fa
samtools dict -o Reference/Ecoli_Ref.dict Reference/Ecoli_Ref.fa
