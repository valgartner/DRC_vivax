#!/usr/bin/env bash
#SBATCH --job-name=biallelic
#SBATCH --mail-type=END
#SBATCH --mail-user=valerie.gartner@duke.edu

module load bcftools
module load htslib/1.3.1-gcb01
module load tabix

VCF=$1

bcftools view -m2 -M2 -v snps ${VCF} > ${VCF%%.vcf.gz}_biallelic_snps_only.vcf 

bgzip ${VCF%%.vcf.gz}_biallelic_snps_only.vcf
tabix -p vcf ${VCF%%.vcf.gz}_biallelic_snps_only.vcf.gz
