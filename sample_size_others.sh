#!/bin/bash
#
#SBATCH --job-name=R_ratio_sample
#SBATCH --output=R_ratio_sample.txt
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --time=1:30:00
#SBATCH -p short-40core

module load shared
module load R/3.6.0

cd /gpfs/scratch/nvitek/molar_ratio_sampling

R

source("molar_ratio_sample_others.R")

