#!/bin/bash
#SBATCH --job-name=bedgraph

# pull out only the chr  & region of interest
# for i in *RR*.dedup_sorted.persitedepth.bedgraph ; do sbatch bedgraph-subregions.sh ${i} ; done

BG=$1


# RBP2a
awk '{if($2>=93324 && $2<=140966) print $0}' ${BG} > ${BG%.dedup*}_RBP2a_10kb-each-side.persitedepth.bedgraph
