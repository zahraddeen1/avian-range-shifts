#!/bin/bash

#SBATCH --job-name=br_lc
#SBATCH --output=br_lc_%A_%a.out
#SBATCH --error=br_lc_%A_%a.err
#SBATCH --array=1-8
#SBATCH --time=6-0
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --mem=10G

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

# Load modules
module load r/4.0.1

# Run rscript
Rscript get_breedingrange_landcover.R
