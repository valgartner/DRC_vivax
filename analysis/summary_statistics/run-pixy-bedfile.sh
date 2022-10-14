#!/bin/bash
#SBATCH --job-name=pixy
#SBATCH --mail-type=END
#SBATCH -c40

source /path/to/miniconda3/bin/activate /path/to/miniconda3/envs/vivax
module load tabix

vcf=$1
popfile=$2

pixy --stats pi dxy --vcf ${vcf} --n_cores 40 --populations ${popfile} --bed-file PVP01.chroms_for-pixy_1kb-windows.bed --output_prefix ${popfile%%.txt}_allchroms_1kb-windows_pixy

