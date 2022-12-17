#!/usr/bin/env Rscript


library(tidyverse)
args = commandArgs(trailingOnly=TRUE) #for script version only

filename <- args[1]

pi_file <- read.table(filename, header=TRUE, skip=0, sep='\t')

grouped_pi_file <- pi_file %>% 
  group_by(pop) %>% 
  drop_na() %>%
  summarise(genome_avg_pi = sum(count_diffs)/sum(count_comparisons))

write_tsv(grouped_pi_file, paste(filename, "_genome-average-pi.txt", sep=""), quote="none")