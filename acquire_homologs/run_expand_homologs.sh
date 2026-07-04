#!/bin/bash
#SBATCH --job-name=ndh2_expand
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=acquire_homologs/logs/expand_%j.out
#SBATCH --error=acquire_homologs/logs/expand_%j.err

set -euo pipefail

mkdir -p acquire_homologs/outputs_expanded acquire_homologs/logs

# Optional: export your NCBI API key for higher rate limits (10 req/s vs 3 req/s)
# export NCBI_API_KEY=...

python -u acquire_homologs/expand_homologs.py acquire_homologs/outputs_expanded

echo "Expansion complete at $(date)"
