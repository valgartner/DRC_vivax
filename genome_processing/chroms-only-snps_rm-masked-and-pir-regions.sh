#!/usr/bin/env bash
#SBATCH --job-name=unmasked_chroms

module load htslib/1.3.1-gcb01 
module load tabix
module load vcftools
module load bcftools

NAME=$1
REF=$2
set -e

#remove masked regions
echo ">>> removing masked regions <<<"
vcftools --gzvcf ${NAME} --out ${NAME}_masked-rm --exclude-bed PVP01.genome.mask.sorted.bed --recode --keep-INFO-all
#out file will look like
#ethiopia-haploid-combined-joint-called_masked-rm.recode.vcf


# chr 10 subtelomere & pir regions
vcftools --vcf ${NAME}_masked-rm.recode.vcf --out ${NAME}_masked-rm.recode.chr10-subtel-rm --exclude-bed chr-10-subtelomeres-from-pi.bed --recode --keep-INFO-all

#extract chromosomes only (no contigs)
echo ">>> extract chromosomes only <<<"
bgzip ${NAME}_masked-rm.recode.chr10-subtel-rm.recode.vcf
bcftools index ${NAME}_masked-rm.recode.chr10-subtel-rm.recode.vcf.gz
bcftools view ${NAME}_masked-rm.recode.chr10-subtel-rm.recode.vcf.gz --regions-file PVP01.chroms.bed > ${NAME}_masked-rm.recode.chr10-subtel-rm.recode_chroms_only.vcf

bgzip ${NAME}_masked-rm.recode.chr10-subtel-rm.recode_chroms_only.vcf
tabix -p vcf ${NAME}_masked-rm.recode.chr10-subtel-rm.recode_chroms_only.vcf.gz

### Uncomment this section to extract SNPs from gVCF
#Extract SNPs
#bcftools view -m2 -M2 -v snps ${NAME}_masked-rm.recode.chr10-subtel-rm.recode_chroms_only.vcf > ${NAME}_masked-rm.recode.chr10-subtel-rm.recode_chroms_only.snps.vcf 

#bgzip ${NAME}_masked-rm.recode.chr10-subtel-rm.recode_chroms_only.snps.vcf
#tabix -p vcf ${NAME}_masked-rm.recode.chr10-subtel-rm.recode_chroms_only.snps.vcf.gz
