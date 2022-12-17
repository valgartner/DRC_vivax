#!/bin/bash
#SBATCH --job-name=bedgraph

# pull out only the chr for each region of interest
# for i in *RR*.dedup_sorted.persitedepth.bedgraph ; do sbatch bedgraph-subregions.sh ${i} ; done

BG=$1


# DBPs

#DBP2
# samtools view -bh ${BAM} "LT635612:94013-117429" > DBP2-extended_${BAM%%_HC_ATTACTCG-CCTATCCT_S5_L001.dedup.bam}.dedup.bam
awk '{if($1=="LT635612") print $0}' ${BG} > ${BG%.bedgraph}_LT635612.persitedepth.bedgraph

#DBP
#samtools view -bh ${BAM} "LT635617:972025-995813" > DBP-extended_${BAM%%_HC_ATTACTCG-CCTATCCT_S5_L001.dedup.bam}.dedup.bam
awk '{if($1=="LT635617") print $0}' ${BG} > ${BG%.bedgraph}_LT635617.persitedepth.bedgraph


# RBPs

# RBP1a, RBP1b
#samtools view -bh ${BAM} "LT635618:61450-90153" > RBP1a-extended_${BAM%%_HC_ATTACTCG-CCTATCCT_S5_L001.dedup.bam}.dedup.bam
awk '{if($1=="LT635618") print $0}' ${BG} > ${BG%.bedgraph}_LT635618.persitedepth.bedgraph


# RBP2a
awk '{if($1=="LT635625") print $0}' ${BG} > ${BG%.bedgraph}_LT635625.persitedepth.bedgraph


# RBP2b
awk '{if($1=="LT635619") print $0}' ${BG} > ${BG%.bedgraph}_LT635619.persitedepth.bedgraph

# RBP2c
awk '{if($1=="LT635616") print $0}' ${BG} > ${BG%.bedgraph}_LT635616.persitedepth.bedgraph


# Other Loci of interest?

# MSP3 gene family
# Ford et al 2020
# PVP01_1031200, PVP01_1031300, PVP01_101400 show signs of selection
# Rice et al 2014

awk '{if($1=="LT635621") print $0}' ${BG} > ${BG%.bedgraph}_LT635621.persitedepth.bedgraph
