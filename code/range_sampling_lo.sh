#!/bin/bash

#SBATCH --job-name=spp_range_lo
#SBATCH --time=5-0
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --mem-per-cpu=10G

# Load modules
module load r/4.0.1
module load geos/3.8.1
module load gdal/3.2.1
module load proj/7.2.1

# Run rscript
Rscript species_range_estimates_longleaf_lo.R