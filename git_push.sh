#!/bin/bash 
#SBATCH -c 8  # Number of Cores per Task
#SBATCH --mem=12g  # Requested Memory
#SBATCH --partition=cpu
#SBATCH --account=pi_alc_umass_edu
#SBATCH -t 4:00:00  # Job time limit
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

# GitHub username and Personal Access Token (PAT)
USERNAME="spconway"
TOKEN="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"

# URL of the repository (remove https:// from here)
REPO_URL="github.com/sp-conway/bw.git"

# Perform the git push using the embedded username and token
git push https://${USERNAME}:${TOKEN}@${REPO_URL}
