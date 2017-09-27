#!/bin/bash

java -jar GenomeAnalysisTK.jar \
	-T GenotypeGVCFs \
	-R Reference/Ecoli_Ref.fa \
	--variant SNPs/SRR098034.gvcf \
	--variant SNPs/SRR098035.gvcf \
	--variant SNPs/SRR098038.gvcf \
	--variant SNPs/SRR098039.gvcf \
	--variant SNPs/SRR098289.gvcf \
	-o SNPs/combined.vcf
