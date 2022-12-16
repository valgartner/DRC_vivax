This repository contains code and steps used for investigating *Plasmodium vivax* evolutionary history in central Africa using whole genome sequencing data. Results of this work will be available on bioRxiv soon.

# *P. vivax* Sample Data

All samples used in this study can be found in the metadata table: `sample_info/metadata_table.csv`.

Multiplicity of Infection (MOI) was determined using the following steps:

1.  generate two vcfs: one polyclone with a max of 3 clones, and one a gvcf which indicates the level of coverage at each base.  The commands to generate these are in `make-fastas-6.sh`
```
$ octopus -I \${DEDUP_BAM} -R ${REF} -T LT635626 -o api.poly3.vcf.gz --annotations AD -C polyclone  --max-clones 3 --threads 16 --sequence-error-model PCR
$ octopus -I \${DEDUP_BAM} -R ${REF} -T LT635626 -o api.g.vcf.gz --annotations AD --refcall POSITIONAL --threads 16 --sequence-error-model PCR
```
2. `check-accessions.py` scans the genomes2 directory that contains directories/files like ERR12355/api.poly3.vcf.gz .  It produced the simple text file genomes2/mono-0.9.txt  with lines like
```
ERR773745       OK
ERR773746       PolyClonal
ERR773747       NoGVCF
ERR773748       OK
```
The 0.9 indicates that a site is considered homozygous if the major allele frequency is 0.9. Run it like `check-accesions.py --cutoff 0.9 --allowed-het-sites=1 > mono-0.9.txt`

Create the list of accessions with `grep OK mono-0.9.txt | cut -f1 > mono-0.9-accs.txt`

Accessions with no MOI data were included by default.

# Variant Calling

All scripts needed to download, map to PvP01 reference genome, and do variant calling are included in the `genome_processing` directory. 

## Download fastq and map to PvP01 reference genome
Run `get-haploid-gvcf.sh` using the command
```
$ sbatch get-haploid-gvcf.sh accessions.txt
``` 
where "accessions.txt" contains one run accession number per line. This script will launch individual bash scripts via an array for each accession number to run in parallel, limited by the array size indicated (currently set to run 10 accessions at one time: `#SBATCH --array=1-${NACC}%10`).

This script will create a new directory for each accession that contains the mapped and deduplicated BAM file as well as the gVCF file. Before combining individual VCF files for the joint callling step, the gVCF file needs to be updated to include the sample name for the genotype information. Run `add-sample-name.sh` with each accession number as the input to update the gVCF file.
```
$ for i in $(< accessions.txt) ; do sbatch add-sample-name.sh ${i} ; done
```
This script will create a new file called `<accesion>-samp.g.vcf` in the same directory as the original gVCF file.

## Combine individual gVCF files
To combine individual gVCFs into one file, first create a list of file locations:
```
$ for i in $(< accessions.xt) ; do echo ${i}/${i}-samp.g.vcf ; done > accessions.list
```
Then run
```
$ sbatch combine-gvcfs.sh accessions.list
```
This script will create a new gVCF file named `accessions-combined.g.vcf.gz`. 

## Joint Calling
Run the joint calling step with the command:
```
$ sbatch joint_call.sh accessions-combined.g.vcf.gz
```
The output of this script will be a file called `accessions-combined-joint-called.g.vcf.gz`

## Remove contigs and hypervariable regions
Before using this gVCF file for analyses, remove unmapped contig reads and reads from hypervariable sites:
keep only chromosomes (no contigs) and remove masked regions using the command:
```
$ sbatch chroms-only-snps_rm-masked-and-pir-regions.sh accessions-combined-joint-called.g.vcf.gz
```
Note: the script `chroms-only-snps_rm-masked-and-pir-regions.sh` can also be used to extract biallelic SNPs by uncommenting the final lines. 



# Analyses

## Figure 1: Admixture Analysis
All scripts can be found in `analysis/admixture_analysis/`.

Starting with SNPs-only VCF file (`/genome_processing/chroms-only-snps_rm-masked-and-pir-regions.sh` set to output biallelic SNPs only):

### Downsample loci
```
$ sort -R min_filt_no-singletons.recode.pruned.genotypes.bim | head -n 100000 | awk '{print $2}' > random100k.snps
#https://www.biostars.org/p/16038/#16085

# change format
# old: chr:pos
# now: chr	pos
$ sed "s/\:/\t/g" random100k.snps > random100k.snps.txt

# extract random positions
$ bcftools view -R random100k.snps.txt min_filt_no-singletons.recode.vcf.gz > random100k_min_filt_no-singletons.recode.vcf

# SORT positions with vcftools 'vcf-sort' tool
$ cat random100k_min_filt_no-singletons.recode.vcf | vcf-sort > random100k-SORTED_min_filt_no-singletons.recode.vcf
```

Convert PvP01 chromosome names to integers
```
$ for i in $(<replace-chr-w-ints_sed-arguments.txt ) ; do sed -i ${i} chr-as-int_global_vivax.vcf ; done
```

### Running AdmixturePipeline

AdmixturePipeline can be downloaded from GitHub: [https://github.com/stevemussmann/admixturePipeline](https://github.com/stevemussmann/admixturePipeline)

To use on haploid data, update the command string in script [admixture.py](https://github.com/stevemussmann/admixturePipeline/blob/master/admixture.py) to indicate haploid mode:
```
#### HAPLOID MODE
command_string = "admixture" + " -s " + str(np.random.randint(1000000)) + " --cv=" + str(self.cv) + " " + self.prefix + ".ped " + str(i) + " " + haploid_str
self.run_program(command_string,i,j)
```

To run, you will need a tab-separated file named `popmap.txt` that contains:
```
accession1      Population1
accession2      Population2
... etc.
```

Submit the pipeline to Slurm using the script:
```
$ sbatch run-admixturePipeline.sh chr-as-int_global_vivax.vcf
```

### Visualize Admixture output using Pong
Pong software can be downloaded from the Github repo: [https://github.com/ramachandran-lab/pong](https://github.com/ramachandran-lab/pong)


Prepare information files for generating Admixture visualization
```
# File map
$ for i in *.Q ; do j=$i ; i=${i##*.pruned.genotypes.} ; echo -e "k-${i%%.Q}\t${i%%_[0-9]*.Q}\t${j}" ; done > filemap.txt

# Pop order
$ awk '{print $2}' popmap.txt |sort|uniq > pop_order.txt
# then manually arrange the order in this file to represent left to right on map	

# Mapping individuals to populations
$ awk '{print $2}' popmap.txt > ind2pop.txt

#Run Pong
$ pong -m filemap.txt -i ind2pop.txt -n pop_order_revised.txt -v
```
Then view results in browser.

To generate Cross Validation Error box plots:
```
for i in *.stdout ; do grep -h "CV error" ${i} >> overall_cv_summary.txt ; done
awk '{print $3"\t"$4}' overall_cv_summary.txt | tr -d '():K=' >> cv_summary_table.txt
```
and visualize with `analysis/admixture_analysis/CrossValidationError_boxplots.Rmd`


## Figure 1: PCA
This analysis uses Plink files generated by AdmixturePipeline in the previous step.
```
$ plink --bfile global_population.genotypes --pca 
```
And visualize with `analysis/plink_pca/global-vivax-pca.Rmd`


## Figure 1: Phylogenetic trees
Bash scripts to start Slurm jobs are included in `analysis/trees`. 

Starting with gvcf file produced by `/genome_processing/chroms-only-snps_rm-masked-and-pir-regions.sh`, first convert VCF to phylip format using script: [vcf2phylip](https://github.com/joanam/scripts/blob/master/vcf2phylip.py)
```
$ sbatch vcf_to_phy.sh vivax.vcf.gz
```

Then remove invariant sites using script from [https://github.com/btmartin721/raxml_ascbias](https://github.com/btmartin721/raxml_ascbias)
```
$ sbatch rm_invariants.sh vivax.phy
```
which will produce a new phylip file `vivax.invariants-rm_snps-only.phy`

To run IQtree, use the command:
```
$ sbatch run_iqtree.sh vivax.invariants-rm_snps-only.phy
```
And visualize using FigTree or software of your choice.

# Supplemental Information

## Figure S1: *P. vivax* genome private alleles as a measure of population variation, separated by continent.

## Figure S2: Genome-wide Nucleotide Diversity within Africa

## Table S1: Summary Statistics
Scripts for this section can be found in `analysis/summary_statistics`

### Segregating Sites
First filter the VCF in two ways:
1. Keep only biallelic sites (sites where every individual either has the reference allele or has a single alternate allele)
```
$ bcftools view -m2 -M2 -v snps ${VCF} > ${VCF%%.vcf.gz}_biallelic_snps_only.vcf 
```

2. Keep only sites where at least one individual in the VCF file has the alternate allele in the GT field (1)
```
$ bcftools view -S ${SUBSAMPLE} ${VCF} --min-ac=1 > ${SUBSAMPLE%%-accessions.txt}_${VCF%%.vcf.gz}_min-ac1.vcf
```

Then calculate the number of segregating sites in the population:
```
$ for i in *_biallelic_snps_only_min-ac1.vcf.gz ; do ./segregating-sites-from-vcf.py --vcf ${i} ; done
```

### Nucleotide Diversity (Pi)

Create window bedfile

`bedtools makewindows -b PVP01.chroms.bed -w 1000 > PVP01.chroms_1kb-windows.bed`

For each population:

`./run-pixy-bedfile.sh`

### Private Alleles
Visualize using `visualize-private-alleles.Rmd`

## Supplementary Table 2: Identification of potential gene duplications in DRC *P. vivax* using read depth.


## Supplemental Figure 4B. Duplication of PvDBP in African samples