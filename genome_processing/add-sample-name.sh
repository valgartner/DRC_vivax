#!/bin/bash
#SBATCH --job-name=addsampname

# to run, provide accession number (and update hard-coded path to gvcf file if necessary)
# e.g. for i in $(< africa-accessions) ; do sbatch add-sample-name.sh ${i} ; done


SAMP=$1
IN=./${SAMP}/${SAMP}.g.vcf
OUT=./${SAMP}/${SAMP}-samp.g.vcf

module load jdk/1.8.0_45-fasrc01 # load java module


RUN_PICARD="java -jar -Xmx7g picard.jar"

${RUN_PICARD} RenameSampleInVcf INPUT=${IN} OUTPUT=${OUT} NEW_SAMPLE_NAME=${SAMP}
