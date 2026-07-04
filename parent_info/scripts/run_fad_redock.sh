#!/bin/bash
#SBATCH --job-name=fad_redock
# partition omitted to use default CPU partition
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/fad_redock_%j.out
#SBATCH --error=logs/fad_redock_%j.err

set -euo pipefail

# Paths
BASE_DIR="/resnick/groups/shapirolab/trao2/ndh2_design/parent_info"
INPUT_PDB="${BASE_DIR}/holo_ndh2_s_aureus_insp.pdb"
OUTDIR="${BASE_DIR}/outputs/fad_refine_${SLURM_JOB_ID:-local}"
export INPUT_PDB
export OUTDIR

mkdir -p "${BASE_DIR}/logs"
mkdir -p "${OUTDIR}"

echo "=== FAD redocking / local refinement ==="
echo "Start: $(date)"
echo "Input: ${INPUT_PDB}"
echo "Out:   ${OUTDIR}"

# Activate docking environment (vina + openbabel expected here)
source /resnick/groups/shapirolab/trao2/miniconda3/bin/activate docking

# 1) Split receptor and ligand (detect FAD by residue name)
python - << 'PY'
import os
pdb = os.environ['INPUT_PDB']
out = os.environ['OUTDIR']
rec_path = os.path.join(out, 'receptor.pdb')
lig_path = os.path.join(out, 'fad.pdb')
with open(pdb) as f, open(rec_path, 'w') as rec, open(lig_path, 'w') as lig:
    for line in f:
        if line.startswith('ATOM'):
            rec.write(line)
        elif line.startswith('HETATM') and line[17:20].strip() == 'FAD':
            lig.write(line)
print(f"Wrote {rec_path} and {lig_path}")
PY

# 2) Prepare PDBQT (protein rigid, ligand with hydrogens)
obabel "${OUTDIR}/receptor.pdb" -O "${OUTDIR}/receptor.pdbqt" -xr
obabel "${OUTDIR}/fad.pdb"      -O "${OUTDIR}/fad.pdbqt" -xh

# 3) Compute docking box centered on original FAD pose
python - << 'PY'
import os, numpy as np
out = os.environ['OUTDIR']
xyz = []
with open(os.path.join(out, 'fad.pdb')) as f:
    for line in f:
        if line.startswith(('ATOM','HETATM')):
            try:
                xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except: pass
arr = np.array(xyz)
minc, maxc = arr.min(0), arr.max(0)
center = (minc + maxc)/2.0
size = (maxc - minc) + 8.0  # 8A margin around original FAD envelope
with open(os.path.join(out,'box.txt'),'w') as fp:
    fp.write(f"center {center[0]} {center[1]} {center[2]}\n")
    fp.write(f"size {size[0]} {size[1]} {size[2]}\n")
print('Box:', center, size)
PY

# 4) Redock FAD with AutoDock Vina (local optimization within box)
python - << 'PY'
import os
from vina import Vina
out = os.environ['OUTDIR']
with open(os.path.join(out,'box.txt')) as f:
    lines = f.read().strip().splitlines()
center = list(map(float, lines[0].split()[1:]))
size   = list(map(float, lines[1].split()[1:]))

v = Vina(sf_name='vina', seed=42, verbosity=1)
v.set_receptor(os.path.join(out,'receptor.pdbqt'))
v.set_ligand_from_file(os.path.join(out,'fad.pdbqt'))
try:
    import os as _os
    cpu = int(_os.environ.get('SLURM_CPUS_PER_TASK','16'))
    v.set_cpu(cpu)
except Exception:
    pass
v.compute_vina_maps(center=center, box_size=size)
v.dock(exhaustiveness=64, n_poses=5)
E = v.energies(n_poses=1)[0][0]
open(os.path.join(out,'fad_energy.txt'),'w').write(f"{E:.2f}\n")
v.write_poses(os.path.join(out,'fad_out.pdbqt'), n_poses=1)
print(f"Vina best energy: {E:.2f} kcal/mol")
PY

# 5) Convert docked pose to PDB and merge into a visualization complex
obabel "${OUTDIR}/fad_out.pdbqt" -O "${OUTDIR}/fad_out.pdb"
{
  echo "REMARK Complex: receptor + redocked FAD";
  grep -E '^(ATOM|HETATM|TER|MODEL|ENDMDL)' "${OUTDIR}/receptor.pdb";
  echo "TER";
  grep -E '^(ATOM|HETATM)' "${OUTDIR}/fad_out.pdb";
  echo "END";
} > "${OUTDIR}/complex_receptor_fad_redocked.pdb"

conda deactivate

echo "End: $(date)"
