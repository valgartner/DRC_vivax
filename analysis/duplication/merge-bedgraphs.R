#!/usr/bin/env Rscript

# before running this script, get each sample in its own file in the format
# sample_acc
# read_depth
# for i in *persitedepth.bedgraph ; do echo -e "${i%%_*}" > ${i%%_*}_depth-col.bedgraph ; awk '{print $3}' ${i} >> ${i%%_*}_depth-col.bedgraph ; done

# run like:
# for i in SANRU*_10kb-each-side.persitedepth_named.bedgraph ; do Rscript merge-bedgraphs.R ${i} ; done
# or
# Rscript merge-bedgraphs.R SANRU_DBP_10kb-each-side.persitedepth_named.bedgraph

library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

sanru_file <- args[1]

sanru <- read.table(sanru_file, header=T, skip=0, sep='\t')

full_table <- sanru

files_list <- list.files(pattern="[A-Z]*RR[0-9]*_depth-col.bedgraph")

for (item in files_list) {
  #print(item)
  new_col <- read.table(item, header=T,skip=0,sep="\t")
  if (nrow(full_table) == nrow(new_col)) {
    full_table <- cbind(full_table, new_col)
  }
}

write_tsv(full_table, "merged.bedgraph", quote="none")
