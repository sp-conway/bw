#!/bin/bash 
#SBATCH -c 24  # Number of Cores per Task
#SBATCH --mem=40g  # Requested Memory
#SBATCH --partition=cpu
#SBATCH --account=pi_alc_umass_edu
#SBATCH -t 24:00:00  # Job time limit
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
module load r-rocker-ml-verse/4.4.0+apptainer
Rscript analysis/crit_analysis_order_constraints.R