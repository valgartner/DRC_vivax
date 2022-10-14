#!/usr/bin/env python3

import sys
import argparse
import re
from pathlib import Path
import vcf
import gzip

def read_accessions(filename):
    acc_file = None
    if filename == '-':
        acc_file = sys.stdin
    else:
        acc_file = open(filename,"r")

    accessions = []
    for line in acc_file:
        line = line.rstrip()
        if line:
            accessions += re.split('\s+', line)
    return accessions

def n_het_sites(poly3_filename, args):
    poly3_reader = vcf.Reader(filename=str(poly3_filename))
    n=0
    for record in poly3_reader:
        data = record.samples[0].data
        GT = record.samples[0].data.GT
        if len(GT) > 1:
            HapFreq = data.MAP_HF
            if max(HapFreq) < args.cutoff:
                n+=1
    return n

        
def check_accession(accession, args):
    dir = Path(accession)
    if not dir.exists():
        return "NoDir"
    gvcf_path = dir / "api.g.vcf.gz"
    if not gvcf_path.exists():
        return "NoGVCF"
    poly3_path = dir / "api.poly3.vcf.gz"
    if not poly3_path.exists():
        return "NoPoly3VCF"

    if n_het_sites(poly3_path, args) > args.allowed_het_sites:
        return "PolyClonal"

    return "OK"


parser = argparse.ArgumentParser(description="Check accessions for MOI")
parser.add_argument("accessions",help="File containing list of accessions")
parser.add_argument("--cutoff",type=float,default=0.95,help="Homozygous alleles must have a frequency greater than this.")
parser.add_argument("--allowed_het_sites",type=float,default=1,help="Maximum number of heterozygous sites allowed")
args = parser.parse_args()

accessions = read_accessions(args.accessions)

for accession in accessions:
    reason = check_accession(accession ,args)
    print(f"{accession}\t{reason}")
