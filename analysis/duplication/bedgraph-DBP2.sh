#!/bin/bash
#SBATCH --job-name=bedgraph

# pull out only the chr  & region of interest
# for i in *RR*.dedup_sorted.persitedepth.bedgraph ; do sbatch bedgraph-subregions.sh ${i} ; done

BG=$1


#DBP2
awk '{if($2>=93013 && $2<=127429)print $0}' ${BG} > ${BG%.dedup*}_DBP2_10kb-each-side.persitedepth.bedgraph

