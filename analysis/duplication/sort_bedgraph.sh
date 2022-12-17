#!/bin/bash
#SBATCH --job-name=sort


module load bedtools2

BED=$1

sort -k 1,1 -k2,2n ${BED} > ${BED%%.bedgraph}.sorted.bedgraph
