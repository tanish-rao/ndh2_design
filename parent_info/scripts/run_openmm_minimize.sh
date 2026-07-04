#!/bin/bash
#SBATCH --job-name=openmm_min_parent_fad
# partition omitted to use default CPU queue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/openmm_min_%j.out
#SBATCH --error=logs/openmm_min_%j.err
#SBATCH --account=shapirolab

set -euo pipefail
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

BASE_DIR="/resnick/groups/shapirolab/trao2/ndh2_design/parent_info"
WORKDIR="${BASE_DIR}/outputs/fad_refine_64538480"

mkdir -p "${BASE_DIR}/logs"

# Initialize conda
source /resnick/groups/shapirolab/trao2/miniconda3/etc/profile.d/conda.sh || true

# Create env if missing
if ! conda env list | grep -q "^openmm_min\b"; then
  conda create -y -n openmm_min python=3.10
fi
conda activate openmm_min

# Install dependencies (idempotent)
conda install -y -c conda-forge openmm openmmforcefields openff-toolkit rdkit openbabel ambertools

# GAFF2 parametrization with AmberTools; first prepare ligand/mol2 with charges
if command -v obabel >/dev/null 2>&1; then
  pushd "${WORKDIR}" >/dev/null
  # Prefer PDBQT (has bond orders from Vina); else PDB
  SRC=""
  if [ -f fad_out.pdbqt ]; then SRC="fad_out.pdbqt"; fi
  if [ -z "$SRC" ] && [ -f fad_out.pdb ]; then SRC="fad_out.pdb"; fi
  if [ -n "$SRC" ]; then
    # Convert to MOL2 for antechamber input
    if [[ "$SRC" == *.pdbqt ]]; then
      obabel -ipdbqt "$SRC" -omol2 -O fad_gaff_input.mol2 || true
    else
      obabel -ipdb   "$SRC" -omol2 -O fad_gaff_input.mol2 || true
    fi
  fi
  popd >/dev/null
fi

# Run antechamber/parmchk2 to create GAFF2 parameters for FAD
pushd "${WORKDIR}" >/dev/null
# Set desired formal charge for FAD (oxidized state)
LIG_CHARGE=-2
# Remove previous GAFF/Amber outputs to ensure regeneration with new charge
rm -f fad_gaff.mol2 fad_gaff.frcmod complex.prmtop complex.inpcrd complex_amber.pdb leap.in leap.log sqm.in sqm.out ANTECHAMBER_* ATOMTYPE.INF >/dev/null 2>&1 || true
if [ -f fad_gaff_input.mol2 ]; then
  echo "[INFO] Running antechamber (AM1-BCC) with charge ${LIG_CHARGE}..."
  if ! antechamber -i fad_gaff_input.mol2 -fi mol2 -o fad_gaff.mol2 -fo mol2 -c bcc -s 2 -nc ${LIG_CHARGE} -at gaff2; then
    echo "[WARN] antechamber AM1-BCC failed; trying 'gas' charges..."
    rm -f fad_gaff.mol2 >/dev/null 2>&1 || true
    if ! antechamber -i fad_gaff_input.mol2 -fi mol2 -o fad_gaff.mol2 -fo mol2 -c gas -s 2 -nc ${LIG_CHARGE} -at gaff2; then
      echo "[WARN] antechamber 'gas' charges failed; trying 'mul' charges..."
      rm -f fad_gaff.mol2 >/dev/null 2>&1 || true
      antechamber -i fad_gaff_input.mol2 -fi mol2 -o fad_gaff.mol2 -fo mol2 -c mul -s 2 -nc ${LIG_CHARGE} -at gaff2 || true
    fi
  fi
  if [ -f fad_gaff.mol2 ]; then
    parmchk2 -i fad_gaff.mol2 -f mol2 -o fad_gaff.frcmod || true
  else
    echo "[ERROR] Failed to create fad_gaff.mol2 with antechamber." >&2
  fi
fi

# Build Amber complex with tleap if ligand params exist
if [ -f fad_gaff.mol2 ] && [ -f fad_gaff.frcmod ] && [ -f receptor.pdb ]; then
  cat > leap.in <<'LEAP'
source leaprc.protein.ff14SB
source leaprc.gaff2
loadamberparams fad_gaff.frcmod
FAD = loadmol2 fad_gaff.mol2
REC = loadpdb receptor.pdb
complex = combine {REC FAD}
set default PBradii mbondi2
saveamberparm complex complex.prmtop complex.inpcrd
savepdb complex complex_amber.pdb
quit
LEAP
  tleap -f leap.in || true
fi
popd >/dev/null
 
# Run minimization
python "${BASE_DIR}/minimize_fad_complex.py" --workdir "${WORKDIR}"

echo "Done at $(date)"
