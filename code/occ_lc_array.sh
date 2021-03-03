#!/bin/bash

#SBATCH --job-name=occ_lc
#SBATCH --output=occ_lc_%A_%a.out
#SBATCH --error=occ_lc_%A_%a.err
#SBATCH --array=1-6
#SBATCH --time=6-0
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --mem=4G

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# Load modules
module load r/4.0.1

# Run rscript
Rscript get_occ_landcover.R
