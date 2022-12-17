#!/usr/bin/env Rscript
# run like
# for i in */ ; do Rscript ./cnv_coords.R ${i%/} ; done 
# where there is one directory for each gene name and all the merged bedgraphs are in the current working dir

library(tidyverse)

args = commandArgs(trailingOnly=TRUE) #for script version only
gene_name <- args[1]

print("XXXXXXXXXXXXXXXXXXXXXXXXXX")
print(gene_name)
print("START")
##### 1. Read in info files
popmap <- read.table("africa_popmap.txt", header=F, skip=0, sep="\t")
popmap <- dplyr::rename(popmap, Accession = V1, Country = V2, Region = V3)

coordinates <- read.table("genes_of_interest_coords.txt", header=T, quote="#", sep="\t")


# select gene of interest for this script
coordinates <- coordinates %>% dplyr::filter(Name == gene_name)

# add info for where the data starts and ends relative to the gene location (start and end refers to gene coordinates)
coordinates <- coordinates %>% 
  mutate(depth_start = start - 10000) %>%
  mutate(depth_end = end + 10000)


##### 2. Read in coverage info
gene_coordinates <- coordinates %>%
  dplyr::filter(Name == gene_name)

bedgraph <- read.table(paste(gene_name, "_merged.bedgraph", sep=""), header=T, skip=0, sep="\t")

##### 3. find coverage levels
# print coordinates to be used
print(paste(gene_coordinates$depth_start, "-", gene_coordinates$start, sep=''))
print(paste(gene_coordinates$start, "-", gene_coordinates$end, sep=''))
print(paste(gene_coordinates$end, "-", gene_coordinates$depth_end, sep=''))


print("END")
print(gene_name)
print("XXXXXXXXXXXXXXXXXXXXXXXXXX")