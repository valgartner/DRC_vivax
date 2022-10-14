#!/bin/bash
#SBATCH --job-name=bgzip
#SBATCH --mail-type=END

# run like: $ ./to-bgzip.sh sample.vcf
# bgzips and indexes the vcf file

module load htslib/1.3.1-gcb01
module load tabix

FILE=$1

bgzip ${FILE}
tabix -p vcf ${FILE}.gz
