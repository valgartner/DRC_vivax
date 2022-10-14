#!/usr/bin/env bash
#SBATCH --mail-type=END
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16

module load python/3.7.4-gcb01
module load plink #pink1.9
module load vcftools
module load #bcftools

infile=$1

admixturePipeline.py -m popmap.txt -v ${infile} -k 2 -K 17 -n 16 -R 50

#k smallest # of populations you want to model
#K biggest # of populations you want to model
#n number of processors for Admixture to use
#a min allele freq
#t Filter loci by thinning out any loci falling within the specified proximity to one another, measured in basepairs. 
#  (default = off, specify an integer greater than 0 to turn it on).
#  Previously tried -t 100
#R number of independant runs for each value of K
#c Specify the cross-validation number for the admixture program. See the admixture program manual for more information (default = 20)
