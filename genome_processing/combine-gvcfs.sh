#!/bin/bash
# Run like `sbatch ./combine-gvcfs.sh southamerica.list`
# NOTE: input file must have the extension '.list' and on each line you should have [ACCESSION].g.vcf to point to the file

#SBATCH --mem=60G
#SBATCH --mail-type=END


module load bwa
module load samtools
module load jdk/1.8.0_45-fasrc01
module load htslib/1.3.1-gcb01
module load tabix 

# Exit immediately if any command returns a failing exit status
set -e

#PvP01 vivax reference genome
reference=PVP01.fasta
ref=${reference}

input=$1

RUN_GATK3="java -jar GenomeAnalysisTK.jar"

function CombineGVCFs ()
{
${RUN_GATK3} -T CombineGVCFs -R ${ref} --variant ${input} -o ${input}-combined.g.vcf
}
CombineGVCFs #call the function

bgzip ${input}-combined.g.vcf
tabix -p vcf ${input}-combined.g.vcf.gz