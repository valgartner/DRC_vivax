#!/usr/bin/env bash
#SBATCH --job-name=bcft_view
#SBATCH --mail-type=END

module load bcftools

VCF=$1
SUBSAMPLE=$2

bcftools view -S ${SUBSAMPLE} ${VCF} > ${SUBSAMPLE%%.txt}.vcf
