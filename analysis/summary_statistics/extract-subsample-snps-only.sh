#!/usr/bin/env bash
#SBATCH --job-name=bcft_view
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu

module load bcftools
module load htslib/1.3.1-gcb01
module load tabix

VCF=$1
SUBSAMPLE=$2

bcftools view -S ${SUBSAMPLE} ${VCF} --min-ac=1 > ${SUBSAMPLE%%-accessions.txt}_${VCF%%.vcf.gz}_min-ac1.vcf
#https://www.biostars.org/p/203809/#204295

bgzip ${SUBSAMPLE%%-accessions.txt}_${VCF%%.vcf.gz}_min-ac1.vcf
tabix -p vcf ${SUBSAMPLE%%-accessions.txt}_${VCF%%.vcf.gz}_min-ac1.vcf.gz
