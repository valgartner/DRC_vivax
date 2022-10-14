#!/usr/bin/env bash
#SBATCH --job-name=bcftools-private-alleles
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu

module load bcftools
module load htslib/1.3.1-gcb01
module load tabix

VCF=$1
ACC_LIST=$2



bcftools view -x --samples-file ${ACC_LIST} ${VCF} > ${ACC_LIST%%-accessions.txt}-private_alleles_biallelic-only.vcf

bgzip ${ACC_LIST%%-accessions.txt}-private_alleles_biallelic-only.vcf
tabix -p vcf ${ACC_LIST%%-accessions.txt}-private_alleles_biallelic-only.vcf.gz
