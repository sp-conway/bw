#!/bin/bash 
#SBATCH -c 60  # Number of Cores per Task
#SBATCH --mem=400g  # Requested Memory
#SBATCH --partition=cpu
#SBATCH --account=pi_alc_umass_edu
#SBATCH -t 48:00:00  # Job time limit
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
module load r-rocker-ml-verse/4.4.0+apptainer
Rscript analysis/bayes_maxdiff.R 
Rscript analysis/bayes_maxdiff_save_p_rank.R
Rscript analysis/bayes_maxdiff_plot_model_means_v_data_means.R
Rscript analysis/bayes_maxdiff_plot_subject_means_v_data_means.R
