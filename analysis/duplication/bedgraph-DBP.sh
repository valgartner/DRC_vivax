#!/bin/bash
#SBATCH --job-name=bedgraph

# pull out only the chr  & region of interest
# for i in *RR*.dedup_sorted.persitedepth.bedgraph ; do sbatch bedgraph-subregions.sh ${i} ; done

BG=$1


#DBP
awk '{if($2>=962025 && $2<=1105813)print $0}' ${BG} > ${BG%.dedup*}_DBP_10kb-each-side.persitedepth.bedgraph
