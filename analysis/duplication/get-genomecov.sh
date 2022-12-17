#!/bin/bash
#SBATCH --job-name=genome_cov


module load bedtools2

BAM=$1

bedtools genomecov -d -ibam ${BAM} > ${BAM%%.bam}.persitedepth.bedgraph

