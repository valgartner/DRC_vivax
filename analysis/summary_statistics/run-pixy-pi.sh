#!/bin/bash
#SBATCH --job-name=pixy
#SBATCH -c40

source /path/to/miniconda3/bin/activate /path/to/miniconda3/envs/vivax
module load tabix

vcf=$1
popfile=$2

pixy --stats pi --vcf ${vcf} --n_cores 40 --populations ${popfile} --window_size 1000 --chromosomes 'LT635612,LT635613,LT635614,LT635615,LT635616,LT635617,LT635618,LT635619,LT635620,LT635621,LT635622,LT635623,LT635624,LT635625,LT635626,LT635627' --output_prefix ${popfile%%.txt}_allchroms_1kb-windows_pixy