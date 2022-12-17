#!/bin/bash
#SBATCH --job-name=bedgraph

# pull out only the chr  & region of interest
# for i in *RR*.dedup_sorted.persitedepth.bedgraph ; do sbatch bedgraph-subregions.sh ${i} ; done

BG=$1

# MSP3 gene family
# Ford et al 2020
# PVP01_1031200, PVP01_1031300, PVP01_101400 show signs of selection
# Rice et al 2014

# MSP3.9
awk '{if($2>=1320308 && $2<=1363631)print $0}' ${BG} > ${BG%.dedup*}_MSP3.9_10kb-each-side.persitedepth.bedgraph    #1340308        1343631 ID=PVP01_1031200;Name=MSP3.9

# MSP3.11
awk '{if($2>=1309792 && $2<=1352916)print $0}' ${BG} > ${BG%.dedup*}_MSP3.11_10kb-each-side.persitedepth.bedgraph  #1329792        1332916 ID=PVP01_1030900;Name=MSP3.11

# MSP3.10
awk '{if($2>=1312979 && $2<=1355525)print $0}' ${BG} > ${BG%.dedup*}_MSP3.10_10kb-each-side.persitedepth.bedgraph  #1332979        1335525 ID=PVP01_1031000;Name=MSP3.10

# MSP3.2
awk '{if($2>=1338155 && $2<=1380860)print $0}' ${BG} > ${BG%.dedup*}_MSP3.2_10kb-each-side.persitedepth.bedgraph   #1358155        1360860 ID=PVP01_1031600;Name=MSP3.2

# MSP3G
awk '{if($2>=1317125 && $2<=1358381)print $0}' ${BG} > ${BG%.dedup*}_MSP3G_10kb-each-side.persitedepth.bedgraph            #1337125        1338381 ID=PVP01_1031100;Name=MSP3G

# MSP3.3
awk '{if($2>=1334475 && $2<=1376829)print $0}' ${BG} > ${BG%.dedup*}_MSP3.3_10kb-each-side.persitedepth.bedgraph   #1354475        1356829 ID=PVP01_1031500;Name=MSP3.3

# MSP3.8
awk '{if($2>=1325019 && $2<=1367715)print $0}' ${BG} > ${BG%.dedup*}_MSP3.8_10kb-each-side.persitedepth.bedgraph   #1345019        1347715 ID=PVP01_1031300;Name=MSP3.8

# MSP3.1
awk '{if($2>=1342254 && $2<=1384791)print $0}' ${BG} > ${BG%.dedup*}_MSP3.1_10kb-each-side.persitedepth.bedgraph   #1362254        1364791 ID=PVP01_1031700;Name=MSP3.1

# MSP3.5
awk '{if($2>=1329145 && $2<=1372525)print $0}' ${BG} > ${BG%.dedup*}_MSP3.5_10kb-each-side.persitedepth.bedgraph   #1349145        1352525 ID=PVP01_1031400;Name=MSP3.5
