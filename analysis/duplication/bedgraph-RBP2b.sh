#!/bin/bash
#SBATCH --job-name=bedgraph

# pull out only the chr  & region of interest
# for i in *RR*.dedup_sorted.persitedepth.bedgraph ; do sbatch bedgraph-subregions.sh ${i} ; done

BG=$1

# RBP2b
awk '{if($2>=13312 && $2<=63704) print $0}' ${BG} > ${BG%.dedup*}_RBP2b_10kb-each-side.persitedepth.bedgraph

#	33312	43704
