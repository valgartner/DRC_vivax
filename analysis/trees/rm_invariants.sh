#!/usr/bin/env bash
#SBATCH --job-name=rm_inv
#SBATCH --mail-type=END
#SBATCH --mem=300G

PHY=$1

source /path/to/miniconda3/bin/activate /path/to/miniconda3/envs/vivax

./ascbias.py -p ${PHY} -o ${PHY%%.phy}.invariants-rm_snps-only.phy
