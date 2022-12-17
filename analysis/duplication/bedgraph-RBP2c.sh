#!/bin/bash
#SBATCH --job-name=bedgraph

# pull out only the chr  & region of interest
# for i in *RR*.dedup_sorted.persitedepth.bedgraph ; do sbatch bedgraph-subregions.sh ${i} ; done

BG=$1


# RBP2c
awk '{if($2>=1438611 && $2<=1487388)print $0}' ${BG} > ${BG%.dedup*}_RBP2c_10kb-each-side.persitedepth.bedgraph

