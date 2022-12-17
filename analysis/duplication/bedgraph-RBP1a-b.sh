#!/bin/bash
#SBATCH --job-name=bedgraph

# pull out only the chr  & region of interest
# for i in *RR*.dedup_sorted.persitedepth.bedgraph ; do sbatch bedgraph-subregions.sh ${i} ; done

BG=$1


# RBP1a
awk '{if(($1="LT635618") && ($2>=60450 && $2<=91153))print $0}' ${BG} > ${BG%.dedup*}_RBP1a_10kb-each-side.persitedepth.bedgraph

# RBP1b
awk '{if(($1="LT635618") && ($2>=49453 && $2<=80211))print $0}' ${BG} > ${BG%.dedup*}_RBP1b_10kb-each-side.persitedepth.bedgraph

