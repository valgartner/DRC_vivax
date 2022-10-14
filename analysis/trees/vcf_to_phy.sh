#!/usr/bin/env bash

#SBATCH --mail-type=END
#SBATCH --job-name=vcf_to_phy
#SBATCH --mem=250G

#source /path/to/miniconda3/bin/activate /path/to/miniconda3/envs/vivax

VCF=$1

./vcf2phylip.py -i ${VCF} -o ${VCF%%.vcf.gz}.phy -m
