#!/bin/bash
#SBATCH --job-name=parent_tmhmm
# partition omitted to use default CPU partition
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=logs/parent_tmhmm_%j.out
#SBATCH --error=logs/parent_tmhmm_%j.err
#SBATCH --account=shapirolab

set -euo pipefail

echo "=== DeepTMHMM for Parent NDH2 ==="
echo "Start: $(date)"

# Ensure logs and output directories exist
cd /resnick/groups/shapirolab/trao2/ndh2_design
mkdir -p logs
mkdir -p parent_info/outputs
# Use a unique output directory (DeepTMHMM requires it to not already exist)
OUTBASE=/resnick/groups/shapirolab/trao2/ndh2_design/parent_info/outputs
OUTDIR=${OUTBASE}/tmhmm_parent_${SLURM_JOB_ID:-$(date +%Y%m%d_%H%M%S)}
echo "Output directory: $OUTDIR"

# Activate DeepTMHMM environment
source /resnick/groups/shapirolab/trao2/miniconda3/bin/activate deeptmhmm_stable

# Run DeepTMHMM prediction
cd /resnick/groups/shapirolab/trao2/tools/deeptmhmm/DeepTMHMM-Academic-License-v1.0
python predict.py \
  --fasta /resnick/groups/shapirolab/trao2/ndh2_design/parent_info/parent_enzyme.fasta \
  --output-dir "$OUTDIR"

echo "End: $(date)"
