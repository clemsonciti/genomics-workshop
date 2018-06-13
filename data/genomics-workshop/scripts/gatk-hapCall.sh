#!/bin/bash
#PBS -N gatk
#PBS -l select=1:ncpus=8:mem=15gb,walltime=2:00:00
#PBS -j oe


echo "START ------------------------------"

module add java/1.8.0 
module load gnu-parallel
module load GATK

src=/home/$USER/genomics-workshop

## Use gnu-parallel to use multiple cores
### within one script
parallel --plus -j 8 java -jar $GATK \
	-T HaplotypeCaller \
	-R $src/Reference/Ecoli_Ref.fa \
	-ploidy 1 \
	-ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 \
	-I {} -o $src/SNPs/{/..}.gvcf ::: $src/BamFiles/*.sort.bam


echo "FINISH ----------------------------"



