#!/usr/bin/env Rscript
# run like
# for i in */ ; do Rscript ./cnv_in_genes.R ${i%/} ; done 
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
print(paste("upstream:", "depth start ", gene_coordinates$depth_start, " depth end ", gene_coordinates$start, sep=''))
print(paste("GENE:", "START ", gene_coordinates$start, " END ", gene_coordinates$end, sep=''))
print(paste("downstream:", "depth start ", gene_coordinates$end, " depth end ", gene_coordinates$depth_end, sep=''))

# non-genic
upstream <- bedgraph %>%
  filter(between(pos, gene_coordinates$depth_start, gene_coordinates$start)) 
downstream <- bedgraph %>%
  filter(between(pos, gene_coordinates$end, gene_coordinates$depth_end))
non_genic_coverage <- rbind(upstream, downstream)
non_genic_mean <- non_genic_coverage %>%
  summarise(across(SANRU:SRR570031, mean))

# genic
genic_coverage <- bedgraph %>%
  filter(between(pos, gene_coordinates$start, gene_coordinates$end))
genic_mean <- genic_coverage %>%
  summarise(across(SANRU:SRR570031, mean))


# summarize ratio
summary <- cbind(round(genic_mean/non_genic_mean, 2)) #2 = round to two digits

# Pivot to long format
non_genic_mean_long <- non_genic_mean %>% pivot_longer(everything(), names_to = "Accession", values_to = "Non-Genic_Coverage")
genic_mean_long <- genic_mean %>% pivot_longer(everything(), names_to = "Accession", values_to = "Genic_Coverage")
summary_long <- summary %>% pivot_longer(everything(), names_to = "Accession", values_to = "Coverage_Ratio")

summary_long <- left_join(summary_long, non_genic_mean_long, by="Accession")
summary_long <- left_join(summary_long, genic_mean_long, by="Accession")

# join with pop file
summary_long <- left_join(summary_long, popmap, by="Accession")

# get grouped averages
summary_country <- summary_long %>% 
  group_by(Country) %>%
  summarize(n= n(), mean_coverage_ratio = round(mean(Coverage_Ratio),2), stdev_coverage_ratio = round(sd(Coverage_Ratio),2))

summary_region <- summary_long %>% 
  group_by(Region) %>%
  summarize(n= n(), mean_coverage_ratio = round(mean(Coverage_Ratio),2), stdev_coverage_ratio = round(sd(Coverage_Ratio),2))

# save summary tables to file
write_tsv(summary_long, paste(gene_name, "_CNV_Africa.txt", sep=""), quote="none")
write_tsv(summary_country, paste(gene_name, "_CNV_Africa_coverage_summary_by_country.txt", sep=""), quote="none")
write_tsv(summary_region, paste(gene_name, "_CNV_Africa_coverage_summary_by_region.txt", sep=""), quote="none")

# visualize
# by country
plot_country <- ggplot(summary_long, aes(x=Country, y=Coverage_Ratio))+ 
  geom_boxplot() +
  geom_jitter() +
  theme_classic() +
  labs(title=paste(gene_name, " CNV in Africa by country", sep=""),x="Population", y = "Genic read depth / non-genic read depth") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_discrete(limits = c("Mauritania", "DRC", "Uganda", "Sudan", "Eritrea", "Ethiopia", "Madagascar"))

ggsave(file=paste(gene_name, "_coverage_by_country.png", sep=""), plot=plot_country, width=12, height=8)

# by region
plot_region <- ggplot(summary_long, aes(x=Region, y=Coverage_Ratio))+ 
  geom_boxplot() +
  geom_jitter() +
  theme_classic() +
  labs(title=paste(gene_name, " CNV in Africa by region", sep=""),x="Population", y = "Genic read depth / non-genic read depth") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_x_discrete(limits = c("East", "Central", "West", "Island"))

ggsave(file=paste(gene_name, "_coverage_by_region.png", sep=""), plot=plot_region, width=12, height=8)

print("END")
print(gene_name)
print("XXXXXXXXXXXXXXXXXXXXXXXXXX")