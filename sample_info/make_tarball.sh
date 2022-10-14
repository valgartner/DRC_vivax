#!/bin/bash
#SBATCH --job-name=tarball
#SBATCH --mail-type=END

# run like:
# $ sbatch ./make_tarball.sh directory

tar_dir=$1

tar -czvf ${tar_dir}.tar.gz ${tar-dir}
