#!/bin/bash
#SBATCH --job-name=ndh2_homologs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=acquire_homologs/logs/fetch_%j.out
#SBATCH --error=acquire_homologs/logs/fetch_%j.err

set -euo pipefail

mkdir -p acquire_homologs/outputs acquire_homologs/logs

# Optional: export your NCBI API key here or rely on environment
# export NCBI_API_KEY=...

# Use system/default Python; adjust env if you have a specific venv
python acquire_homologs/fetch_homologs_ncbi.py acquire_homologs/outputs

echo "Done at $(date)"
