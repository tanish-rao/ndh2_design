"""
Tier 6 Filter: Molecular Docking with NADH and DHNA
Filters sequences by predicted binding affinity to NDH-2 substrates.
- Uses AutoDock Vina or GNINA for docking
- Requires predicted 3D structure (from ESMFold PDB output)
- Substrates: NADH (cofactor) and DHNA (demethylmenaquinone, electron acceptor)
- Filter: binding energy < -7.0 kcal/mol for both substrates
Runs in days on ~200-500 sequences.

Prerequisites:
  - ESMFold PDB structures (from Tier 3)
  - AutoDock Vina or GNINA installed
  - Substrate SDF/PDBQT files for NADH and DHNA
  - AutoDockTools for receptor preparation
"""
import argparse
import subprocess
import os
import json
from pathlib import Path
from Bio import SeqIO


NADH_SMILES = "C1C=CN(C=C1C(=O)N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O"
DHNA_SMILES = "C1=CC=C2C(=C1)C(=CC(=C2O)C(=O)O)O"


def prepare_ligand(smiles, output_pdbqt, ligand_name="LIG"):
    """Convert SMILES to PDBQT using obabel."""
    sdf_file = output_pdbqt.replace(".pdbqt", ".sdf")
    pdb_file = output_pdbqt.replace(".pdbqt", ".pdb")
    
    # Generate 3D structure from SMILES
    result = subprocess.run([
        "obabel", f"-:{smiles}", "-O", sdf_file,
        "--gen3d", "--best", "-h"
    ], capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Warning: obabel 3D generation: {result.stderr}")
    
    # Convert to PDB first
    subprocess.run([
        "obabel", sdf_file, "-O", pdb_file, "-h"
    ], capture_output=True, text=True)
    
    # Convert to PDBQT
    subprocess.run([
        "obabel", pdb_file, "-O", output_pdbqt, "-xh"
    ], capture_output=True, text=True)
    
    return output_pdbqt


def prepare_receptor(pdb_file, output_pdbqt):
    """Prepare receptor PDBQT using obabel."""
    # Convert PDB to PDBQT format
    result = subprocess.run([
        "obabel", pdb_file, "-O", output_pdbqt,
        "-xr"  # Rigid receptor
    ], capture_output=True, text=True)
    
    if result.returncode != 0:
        raise Exception(f"Receptor preparation failed: {result.stderr}")
    
    return output_pdbqt


def run_vina(receptor_pdbqt, ligand_pdbqt, output_pdbqt, center, box_size=(25, 25, 25)):
    """Run AutoDock Vina docking using Python API."""
    from vina import Vina
    
    v = Vina(sf_name='vina')
    
    # Set receptor and ligand
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    
    # Set search space
    cx, cy, cz = center
    sx, sy, sz = box_size
    v.compute_vina_maps(center=[cx, cy, cz], box_size=[sx, sy, sz])
    
    # Dock
    v.dock(exhaustiveness=16, n_poses=5)
    
    # Get best energy (first pose)
    energies = v.energies()
    if len(energies) > 0:
        best_energy = energies[0][0]  # First pose, binding affinity
    else:
        best_energy = 0.0
    
    # Write output
    v.write_poses(output_pdbqt, n_poses=5, overwrite=True)
    
    return best_energy


def get_binding_site_center(pdb_file):
    """Estimate binding site center from protein geometric center."""
    coords = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    coords.append((x, y, z))
                except ValueError:
                    continue
    if not coords:
        return (0.0, 0.0, 0.0)
    cx = sum(c[0] for c in coords) / len(coords)
    cy = sum(c[1] for c in coords) / len(coords)
    cz = sum(c[2] for c in coords) / len(coords)
    return (cx, cy, cz)


def filter_by_docking(input_fasta, output_fasta, structures_dir,
                       min_nadh_energy=-7.0, min_dhna_energy=-7.0,
                       scores_file=None, work_dir="docking_work"):
    """
    Filter sequences by molecular docking scores.
    Expects ESMFold PDB structures in structures_dir/{seq_id}.pdb
    """
    Path(work_dir).mkdir(parents=True, exist_ok=True)
    sequences = {r.id: r for r in SeqIO.parse(input_fasta, "fasta")}

    # Prepare ligands once
    nadh_pdbqt = os.path.join(work_dir, "nadh.pdbqt")
    dhna_pdbqt = os.path.join(work_dir, "dhna.pdbqt")
    print("Preparing ligands...")
    prepare_ligand(NADH_SMILES, nadh_pdbqt, "NADH")
    prepare_ligand(DHNA_SMILES, dhna_pdbqt, "DHNA")

    passed, failed, scores = [], [], {}
    for i, (seq_id, record) in enumerate(sequences.items()):
        pdb_file = os.path.join(structures_dir, f"{seq_id}.pdb")
        if not os.path.exists(pdb_file):
            print(f"  Warning: No structure for {seq_id}, skipping")
            failed.append(record)
            continue

        if i % 10 == 0:
            print(f"  Docking {i+1}/{len(sequences)}: {seq_id}...")

        try:
            receptor_pdbqt = os.path.join(work_dir, f"{seq_id}_receptor.pdbqt")
            prepare_receptor(pdb_file, receptor_pdbqt)
            center = get_binding_site_center(pdb_file)

            nadh_out = os.path.join(work_dir, f"{seq_id}_nadh_out.pdbqt")
            dhna_out = os.path.join(work_dir, f"{seq_id}_dhna_out.pdbqt")

            nadh_energy = run_vina(receptor_pdbqt, nadh_pdbqt, nadh_out, center)
            dhna_energy = run_vina(receptor_pdbqt, dhna_pdbqt, dhna_out, center)

            scores[seq_id] = {"nadh": nadh_energy, "dhna": dhna_energy}

            if nadh_energy <= min_nadh_energy and dhna_energy <= min_dhna_energy:
                passed.append(record)
            else:
                failed.append(record)

        except Exception as e:
            print(f"  Warning: Docking failed for {seq_id}: {e}")
            failed.append(record)

    SeqIO.write(passed, output_fasta, "fasta")

    if scores_file:
        with open(scores_file, "w") as f:
            f.write("sequence_id,nadh_energy,dhna_energy\n")
            for seq_id, s in scores.items():
                f.write(f"{seq_id},{s['nadh']:.2f},{s['dhna']:.2f}\n")
        print(f"  Docking scores saved to {scores_file}")

    total = len(sequences)
    print(f"Molecular docking filter (NADH<{min_nadh_energy}, DHNA<{min_dhna_energy} kcal/mol):")
    print(f"  Input:  {total} sequences")
    print(f"  Passed: {len(passed)} ({100*len(passed)/total:.1f}%)")
    print(f"  Failed: {len(failed)} sequences")
    print(f"  Output: {output_fasta}")
    return passed


def main():
    parser = argparse.ArgumentParser(description="Filter sequences by molecular docking scores.")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    parser.add_argument("--structures_dir", required=True, help="Directory with ESMFold PDB structures")
    parser.add_argument("--min_nadh_energy", type=float, default=-7.0, help="Min NADH binding energy kcal/mol (default: -7.0)")
    parser.add_argument("--min_dhna_energy", type=float, default=-7.0, help="Min DHNA binding energy kcal/mol (default: -7.0)")
    parser.add_argument("--scores_file", default=None, help="Optional: save docking scores to CSV")
    parser.add_argument("--work_dir", default="docking_work", help="Working directory for docking files")
    args = parser.parse_args()

    filter_by_docking(
        args.input, args.output, args.structures_dir,
        args.min_nadh_energy, args.min_dhna_energy,
        args.scores_file, args.work_dir
    )


if __name__ == "__main__":
    main()
