#!/bin/bash
#SBATCH --job-name=parent_docking_fad
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/parent_analysis_fad_%j.out
#SBATCH --error=logs/parent_analysis_fad_%j.err

set -euo pipefail

BASE_DIR="/resnick/groups/shapirolab/trao2/ndh2_design"
RECEPTOR_PDB="${BASE_DIR}/parent_info/outputs/fad_refine_64538480/complex_receptor_fad_minimized.pdb"
OUT_DIR="${BASE_DIR}/pipeline_results/parent_analysis"
mkdir -p "${BASE_DIR}/logs" "${OUT_DIR}"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

# Activate docking environment
source /resnick/groups/shapirolab/trao2/miniconda3/bin/activate docking

# Prepare receptor PDBQT (rigid receptor)
obabel "${RECEPTOR_PDB}" -O "${OUT_DIR}/receptor_fad.pdbqt" -xr

# Prepare ligands (NADH and DHNA)
NADH_SMILES="C1C=CN(C=C1C(=O)N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O"
DHNA_SMILES="C1=CC=C2C(=C1)C(=CC(=C2O)C(=O)O)O"

obabel -:"${NADH_SMILES}" -O "${OUT_DIR}/nadh.pdb" --gen3d
obabel "${OUT_DIR}/nadh.pdb" -O "${OUT_DIR}/nadh.pdbqt" -xh

obabel -:"${DHNA_SMILES}" -O "${OUT_DIR}/dhna.pdb" --gen3d
obabel "${OUT_DIR}/dhna.pdb" -O "${OUT_DIR}/dhna.pdbqt" -xh

# Python docking with Vina API
python -u - <<'PY'
from vina import Vina
import numpy as np
import os
out_dir = os.path.join('/resnick/groups/shapirolab/trao2/ndh2_design','pipeline_results','parent_analysis')
receptor_pdb = os.path.join(out_dir,'receptor_fad.pdbqt')
receptor_pdb_coords = os.path.join(out_dir,'receptor_fad.pdbqt')

# Compute box from receptor PDBQT (fallback to original PDB if needed)
coords=[]
try:
    with open(receptor_pdb) as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    x,y,z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    coords.append([x,y,z])
                except: pass
except FileNotFoundError:
    with open(os.path.join(out_dir,'receptor_fad.pdb')) as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    x,y,z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    coords.append([x,y,z])
                except: pass

coords = np.array(coords)
center = coords.mean(axis=0)
size = coords.max(axis=0) - coords.min(axis=0) + 10

cpu = int(os.environ.get('SLURM_CPUS_PER_TASK','16'))

for lig in ['nadh','dhna']:
    v = Vina(sf_name='vina', seed=42, verbosity=0)
    try:
        v.set_cpu(cpu)
    except Exception:
        pass
    v.set_receptor(receptor_pdb)
    v.set_ligand_from_file(os.path.join(out_dir, f'{lig}.pdbqt'))
    v.compute_vina_maps(center=center.tolist(), box_size=size.tolist())
    v.dock(exhaustiveness=32, n_poses=5)
    energy = v.energies(n_poses=1)[0][0]
    v.write_poses(os.path.join(out_dir, f'parent_{lig}_out.pdbqt'), n_poses=1)
    print(f'{lig.upper()} binding energy: {energy:.2f} kcal/mol')
    with open(os.path.join(out_dir, f'{lig}_energy.txt'), 'w') as fh:
        fh.write(f'{energy:.2f}\n')
PY

echo "Done at $(date)"
