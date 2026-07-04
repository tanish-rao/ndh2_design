#!/bin/bash
#SBATCH --job-name=parent_docking
#SBATCH --partition=sunshine
#SBATCH --gres=gpu:v100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=logs/parent_analysis_%j.out
#SBATCH --error=logs/parent_analysis_%j.err

echo "=========================================="
echo "Parent Enzyme Structure & Docking Analysis"
echo "=========================================="
echo "Start time: $(date)"
echo ""

cd /resnick/groups/shapirolab/trao2/ndh2_design
mkdir -p logs
mkdir -p pipeline_results/parent_analysis
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

# Step 1: Generate ESMFold structure
echo "Step 1: Generating ESMFold structure for parent enzyme..."
source /resnick/groups/shapirolab/trao2/miniconda3/bin/activate my_esmfold
python -u parent_info/analyze_parent_enzyme.py

if [ ! -f "pipeline_results/parent_analysis/parent_ndh2_lactobacillus.pdb" ]; then
    echo "ERROR: Structure generation failed"
    exit 1
fi

echo ""
echo "Step 2: Running molecular docking with NADH and DHNA..."

# Switch to docking environment
conda deactivate
source /resnick/groups/shapirolab/trao2/miniconda3/bin/activate docking

# Prepare receptor (convert PDB to PDBQT)
echo "  Preparing receptor..."
obabel pipeline_results/parent_analysis/parent_ndh2_lactobacillus.pdb \
    -O pipeline_results/parent_analysis/parent_ndh2_lactobacillus.pdbqt \
    -xr

# Prepare ligands (NADH and DHNA)
NADH_SMILES="C1C=CN(C=C1C(=O)N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O"
DHNA_SMILES="C1=CC=C2C(=C1)C(=CC(=C2O)C(=O)O)O"

echo "  Preparing NADH ligand..."
obabel -:"$NADH_SMILES" -O pipeline_results/parent_analysis/nadh.pdb --gen3d
obabel pipeline_results/parent_analysis/nadh.pdb -O pipeline_results/parent_analysis/nadh.pdbqt -xh

echo "  Preparing DHNA ligand..."
obabel -:"$DHNA_SMILES" -O pipeline_results/parent_analysis/dhna.pdb --gen3d
obabel pipeline_results/parent_analysis/dhna.pdb -O pipeline_results/parent_analysis/dhna.pdbqt -xh

# Run docking with AutoDock Vina
echo "  Running NADH docking..."
python -c "
from vina import Vina
import numpy as np
import os
import os

# Load receptor
v = Vina(sf_name='vina', seed=42, verbosity=0)
try:
    cpu = int(os.environ.get('SLURM_CPUS_PER_TASK', '16'))
    v.set_cpu(cpu)
except Exception:
    pass
v.set_receptor('pipeline_results/parent_analysis/parent_ndh2_lactobacillus.pdbqt')

# Calculate search box (entire protein)
coords = []
with open('pipeline_results/parent_analysis/parent_ndh2_lactobacillus.pdb') as f:
    for line in f:
        if line.startswith('ATOM'):
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coords.append([x, y, z])
            except: pass

coords = np.array(coords)
center = coords.mean(axis=0)
size = coords.max(axis=0) - coords.min(axis=0) + 10

v.set_ligand_from_file('pipeline_results/parent_analysis/nadh.pdbqt')
v.compute_vina_maps(center=center.tolist(), box_size=size.tolist())

# Dock
v.dock(exhaustiveness=32, n_poses=5)
nadh_energy = v.energies(n_poses=1)[0][0]

# Save
v.write_poses('pipeline_results/parent_analysis/parent_nadh_out.pdbqt', n_poses=1)

print(f'NADH binding energy: {nadh_energy:.2f} kcal/mol')
with open('pipeline_results/parent_analysis/nadh_energy.txt', 'w') as f:
    f.write(f'{nadh_energy:.2f}\n')
"

echo "  Running DHNA docking..."
python -c "
from vina import Vina
import numpy as np

# Load receptor
v = Vina(sf_name='vina', seed=42, verbosity=0)
try:
    cpu = int(os.environ.get('SLURM_CPUS_PER_TASK', '16'))
    v.set_cpu(cpu)
except Exception:
    pass
v.set_receptor('pipeline_results/parent_analysis/parent_ndh2_lactobacillus.pdbqt')

# Calculate search box
coords = []
with open('pipeline_results/parent_analysis/parent_ndh2_lactobacillus.pdb') as f:
    for line in f:
        if line.startswith('ATOM'):
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coords.append([x, y, z])
            except: pass

coords = np.array(coords)
center = coords.mean(axis=0)
size = coords.max(axis=0) - coords.min(axis=0) + 10

v.set_ligand_from_file('pipeline_results/parent_analysis/dhna.pdbqt')
v.compute_vina_maps(center=center.tolist(), box_size=size.tolist())

# Dock
v.dock(exhaustiveness=32, n_poses=5)
dhna_energy = v.energies(n_poses=1)[0][0]

# Save
v.write_poses('pipeline_results/parent_analysis/parent_dhna_out.pdbqt', n_poses=1)

print(f'DHNA binding energy: {dhna_energy:.2f} kcal/mol')
with open('pipeline_results/parent_analysis/dhna_energy.txt', 'w') as f:
    f.write(f'{dhna_energy:.2f}\n')
"

echo ""
echo "=========================================="
echo "Analysis Complete"
echo "End time: $(date)"
echo "=========================================="

# Display results
if [ -f "pipeline_results/parent_analysis/nadh_energy.txt" ] && [ -f "pipeline_results/parent_analysis/dhna_energy.txt" ]; then
    echo ""
    echo "Parent Enzyme Binding Energies:"
    echo "  NADH: $(cat pipeline_results/parent_analysis/nadh_energy.txt | tr -d '\n') kcal/mol"
    echo "  DHNA: $(cat pipeline_results/parent_analysis/dhna_energy.txt | tr -d '\n') kcal/mol"
    echo ""
    echo "Results saved in: pipeline_results/parent_analysis/"
fi
