#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --mail-type=END

module load htslib/1.3.1-gcb01 
module load tabix
module load vcftools

vcf_1=$1
vcf_2=$2

vcf-merge ${vcf_1} ${vcf_2} > merged.g.vcf

bgzip merged.g.vcf
tabix -p vcf merged.g.vcf.gz
