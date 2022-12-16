#!/bin/bash
#SBATCH --mem=60G
#SBATCH -c16
#SBATCH --job-name=joint-call
#SBATCH --mail-type=END

#run like ./joint-calling.sh population.g.vcf

#https://github.com/broadinstitute/gatk-docs/blob/master/gatk3-tutorials/(howto)_Discover_variants_with_GATK_-_A_GATK_Workshop_Tutorial.md

module load jdk/1.8.0_45-fasrc01
module load htslib/1.3.1-gcb01 
module load tabix

# Exit immediately if any command returns a failing exit status
set -e

reference=PVP01.fa
ref=${reference}

VCF=$1 #the combined gVCF file
OUT=${VCF%%.g.vcf}-joint-called.g.vcf

echo 'VCF file is' ${VCF}
echo 'output file is' ${OUT}


GATK3_JAR=GenomeAnalysisTK.jar
RUN_GATK3="java -jar ${GATK3_JAR}"


${RUN_GATK3} -T GenotypeGVCFs -R ${ref} --variant ${VCF} -allSites -o ${OUT}


bgzip ${OUT}
tabix -p vcf ${OUT}.gz