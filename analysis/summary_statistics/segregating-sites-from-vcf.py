#!/usr/bin/env python3

import numpy as np
import pandas as pd
import gzip
import argparse
import csv
import os.path


"""
This script takes a gzipped VCF file as input and calculates the number of segregating sites in the file. This input VCF file should be filtered in two ways: 

1. Keep only biallelic sites (sites where every individual either has the reference allele or has a single alternate allele)
e.g.
bcftools view -m2 -M2 -v snps ${VCF} > ${VCF%%.vcf.gz}_biallelic_snps_only.vcf 

2. Keep only sites where at least one individual in the VCF file has the alternate allele in the GT field (1)
e.g.
bcftools view -S ${SUBSAMPLE} ${VCF} --min-ac=1 > ${SUBSAMPLE%%-accessions.txt}_${VCF%%.vcf.gz}_min-ac1.vcf

This script will then determine the number of segregating sites by excluding any positions where all individuals in the file have the alternate allele (1).
We exclude such sites because they do not represent variation within the population (only that the population differs from the reference genome at this position).

The output will create or add to a CSV file the population and the segregating sites count.

Run like:
for i in *_biallelic_snps_only_min-ac1.vcf.gz ; do ./segregating-sites-from-vcf.py --vcf ${i} ; done
"""

parser = argparse.ArgumentParser(description="reads in VCF using pandas")
parser.add_argument("--vcf", help="VCF file to import (including gz extension)") # args.vcf
args = parser.parse_args()

infile = args.vcf
pop_name = infile.split("_",1)[0] # get the population name from VCF file input argument

"""
Define functions
"""
### get column names
def get_vcf_col_names(vcf_path): # https://www.biostars.org/p/416324/#9480044
	with gzip.open(vcf_path, "rt") as ifile:
		for line in ifile:
			if line.startswith("#CHROM"):
				vcf_names = [x for x in line.split('\t')]
				  
				#remove newline character
				vcf_names_2 = []
				for name in vcf_names:
					vcf_names_2.append(name.strip())
				break
	ifile.close()
	return vcf_names_2

# 
def s_count(series):
	S = 0
	for row in series:
		if row >= 1:
			S += 1
		else:
			pass
	return(S)

"""
Read in file
"""
### store column names for later use
col_names = get_vcf_col_names(infile)
sample_list = col_names[9:]
sample_count = len(sample_list)

### Read in the file
vcf = pd.read_table(infile, sep="\t", comment='#', header=None, names=col_names)

### split GT column within each sample column
for samp in sample_list:
	vcf[samp] = vcf[samp].str.split(":",1).str[0] # keep only GT information for each column 
	#print(f"sample: {vcf[samp]}")


"""
Count segregating sites.

If the row contains "0" (any individuals with reference allele), add that site to S count
Don't count rows where all samples have GT="1" (alt allele) OR if all samples are a mix of "1" and "." (missing allele)
"""
# drop info columns except for GT info for each sample
vcf.drop(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"], axis=1, inplace=True)

# https://stackoverflow.com/a/29229653/10176950
# all samples in a row are non-zero (e.g. alt or missing) => sum will be 0

if sample_count > 1: 
	rows_with_ref = (vcf == "0").sum(axis=1)

	# If there are any samples in the row with 0 (ref allele), then the sum of boolean operation (vcf == "0").sum(axis=1) will be >= 1. 
	# These are the segregating sites we want to count.
	# We have to leave GT values as strings to accomodate missing GT (".")
	s_to_keep = s_count(rows_with_ref)

# If there is only one sample in the population, keep only sites with alt GT (e.g. "1").
# Don't count missing GTs (e.g. ".")
elif sample_count == 1:	
	rows_with_alt = (vcf == "1").sum(axis=1)

	s_to_keep = s_count(rows_with_alt)
else:
	print(f'check for problems with file {infile}')



"""
Output results to csv
"""

outfilename = "updated_population_seg-sites_count.csv"
file_exists = os.path.isfile(outfilename)
# 'a' means results will be appended to outfilename if it exists
with open(outfilename, 'a') as out_file:
	headers = ['country', 'sample_no', 'segregating_sites']
	out_writer = csv.writer(out_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

	if not file_exists:  # if true, the file doesn't exist yet -> write the header
		out_writer.writerow(headers)

	out_writer.writerow([pop_name, sample_count, s_to_keep])

out_file.close()