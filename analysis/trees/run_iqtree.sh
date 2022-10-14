#!/usr/bin/env bash
#SBATCH --job-name=iqtree
#SBATCH --mail-type=END
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16

# genome-wide snp trees with ascertainment bias correction model

infile=$1
iqtree -s ${infile} -m GTR+ASC -alrt 1000 -bb 1000 -nt AUTO
