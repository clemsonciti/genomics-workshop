#!/bin/bash

java -jar GenomeAnalysisTK.jar -T HaplotypeCaller \
	-I BamFiles/SRR098034.sort.bam \
	-R Reference/Ecoli_Ref.fa \
	-o SNPs/SRR098034.gvcf \
	-ERC GVCF \
	-variant_index_type LINEAR \
	-variant_index_parameter 128000 \
	-ploidy 1
