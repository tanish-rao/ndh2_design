"""
Tier 7 Filter: Density Functional Theory (DFT)
Filters sequences by electronic structure properties of the FAD/NAD binding site.
- Uses Gaussian or ORCA for DFT calculations on binding site residues
- Key metrics: HOMO-LUMO gap, charge transfer, binding site electrostatics
- Filter: HOMO-LUMO gap > threshold, favorable charge distribution
Runs in weeks on ~50-100 sequences (HPC cluster required).

Prerequisites:
  - Gaussian 16 or ORCA installed
  - ESMFold PDB structures (from Tier 3)
  - Docking poses (from Tier 6) to define binding site geometry
  - Python: cclib for parsing DFT output
"""
import argparse
import subprocess
import os
import json
from pathlib import Path
from Bio import SeqIO


# Key residues for NDH-2 FAD/NAD binding site (based on L. plantarum seed)
# These are approximate - refine based on structural alignment
FAD_BINDING_RESIDUES = ["GGxGxxG", "RxY"]  # Rossmann fold motifs
NAD_BINDING_RESIDUES = ["GGxGxxG"]


def extract_binding_site(pdb_file, docking_pose_pdb, radius=6.0):
    """
    Extract residues within radius of docking pose for DFT calculation.
    Returns list of residue coordinates.
    """
    binding_site_atoms = []
    # Load docking pose center
    pose_coords = []
    with open(docking_pose_pdb) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    pose_coords.append((x, y, z))
                except ValueError:
                    continue

    if not pose_coords:
        return []

    cx = sum(c[0] for c in pose_coords) / len(pose_coords)
    cy = sum(c[1] for c in pose_coords) / len(pose_coords)
    cz = sum(c[2] for c in pose_coords) / len(pose_coords)

    # Extract protein atoms within radius of pose center
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    dist = ((x-cx)**2 + (y-cy)**2 + (z-cz)**2)**0.5
                    if dist <= radius:
                        binding_site_atoms.append(line)
                except ValueError:
                    continue

    return binding_site_atoms


def write_orca_input(atoms, output_file, charge=0, multiplicity=None,
                     functional="B3LYP", basis="6-31G*"):
    """Write ORCA DFT input file for binding site."""
    # Count electrons to determine multiplicity
    electron_count = {'C': 6, 'N': 7, 'O': 8, 'S': 16, 'H': 1}
    total_electrons = 0
    
    for atom_line in atoms:
        element = atom_line[76:78].strip() or atom_line[12:16].strip()[0]
        total_electrons += electron_count.get(element, 0)
    
    # Adjust for charge
    total_electrons -= charge
    
    # Set multiplicity based on electron count
    if multiplicity is None:
        if total_electrons % 2 == 0:
            multiplicity = 1  # Singlet for even electrons
        else:
            multiplicity = 2  # Doublet for odd electrons
    
    with open(output_file, "w") as f:
        f.write(f"! {functional} {basis}\n")
        f.write(f"! RIJCOSX def2/J\n\n")
        f.write(f"%maxcore 4000\n")
        f.write(f"%pal nprocs 1 end\n\n")
        f.write(f"* xyz {charge} {multiplicity}\n")
        for atom_line in atoms:
            element = atom_line[76:78].strip() or atom_line[12:16].strip()[0]
            x = float(atom_line[30:38])
            y = float(atom_line[38:46])
            z = float(atom_line[46:54])
            f.write(f"  {element}  {x:.4f}  {y:.4f}  {z:.4f}\n")
        f.write("*\n")


def parse_orca_output(output_file):
    """Parse ORCA output for HOMO-LUMO gap and other properties."""
    import re
    
    homo_energy = None
    lumo_energy = None

    with open(output_file) as f:
        content = f.read()

    # Find the ORBITAL ENERGIES section for SPIN UP (or just ORBITAL ENERGIES for closed shell)
    # Look for pattern: orbital_number  occupancy  energy(Eh)  energy(eV)
    # HOMO is last orbital with OCC=1.0000, LUMO is first with OCC=0.0000
    
    # Find all orbital lines with format: number  occupancy  energy_Eh  energy_eV
    orbital_pattern = r'^\s+(\d+)\s+([\d.]+)\s+([-\d.]+)\s+([-\d.]+)'
    
    orbitals = []
    in_orbital_section = False
    
    for line in content.split('\n'):
        if 'ORBITAL ENERGIES' in line or 'SPIN UP ORBITALS' in line:
            in_orbital_section = True
            continue
        if 'SPIN DOWN ORBITALS' in line:
            break  # Only use spin up for unrestricted
        
        if in_orbital_section:
            match = re.match(orbital_pattern, line)
            if match:
                orb_num = int(match.group(1))
                occupancy = float(match.group(2))
                energy_eV = float(match.group(4))
                orbitals.append((orb_num, occupancy, energy_eV))
    
    # Find HOMO (last occupied) and LUMO (first unoccupied)
    if orbitals:
        occupied = [(n, e) for n, occ, e in orbitals if occ > 0.5]
        unoccupied = [(n, e) for n, occ, e in orbitals if occ < 0.5]
        
        if occupied and unoccupied:
            homo_energy = occupied[-1][1]  # Last occupied
            lumo_energy = unoccupied[0][1]  # First unoccupied
            gap = lumo_energy - homo_energy
            
            return {"homo": homo_energy, "lumo": lumo_energy, "gap": gap}
    
    return None


def filter_by_dft(input_fasta, output_fasta, structures_dir, docking_dir,
                   min_homo_lumo_gap=2.0, scores_file=None, work_dir="dft_work"):
    """
    Filter sequences by DFT-computed electronic properties of binding site.
    min_homo_lumo_gap: minimum HOMO-LUMO gap in eV (larger = more stable)
    """
    Path(work_dir).mkdir(parents=True, exist_ok=True)
    sequences = {r.id: r for r in SeqIO.parse(input_fasta, "fasta")}

    passed, failed, scores = [], [], {}
    for i, (seq_id, record) in enumerate(sequences.items()):
        pdb_file = os.path.join(structures_dir, f"{seq_id}.pdb")
        nadh_pose = os.path.join(docking_dir, f"{seq_id}_nadh_out.pdbqt")

        if not os.path.exists(pdb_file) or not os.path.exists(nadh_pose):
            print(f"  Warning: Missing structure or docking pose for {seq_id}, skipping")
            print(f"    Looking for: {pdb_file}")
            print(f"    Looking for: {nadh_pose}")
            failed.append(record)
            continue

        print(f"  DFT calculation {i+1}/{len(sequences)}: {seq_id}...")

        try:
            # Extract binding site
            binding_atoms = extract_binding_site(pdb_file, nadh_pose)
            if not binding_atoms:
                print(f"  Warning: No binding site atoms found for {seq_id}")
                failed.append(record)
                continue

            # Write ORCA input
            orca_input = os.path.join(work_dir, f"{seq_id}.inp")
            orca_output = os.path.join(work_dir, f"{seq_id}.out")
            write_orca_input(binding_atoms, orca_input)

            # Run ORCA (use full path for parallel runs)
            orca_exe = "/central/software9/external/orca/orca_6_1_0_linux_x86-64_shared_openmpi418_avx2/orca"
            with open(orca_output, "w") as out_f:
                subprocess.run([orca_exe, orca_input], stdout=out_f, stderr=subprocess.STDOUT)

            # Parse results
            dft_results = parse_orca_output(orca_output)
            if dft_results is None:
                print(f"  Warning: DFT parsing failed for {seq_id}")
                failed.append(record)
                continue

            scores[seq_id] = dft_results
            if dft_results["gap"] >= min_homo_lumo_gap:
                passed.append(record)
            else:
                failed.append(record)

        except Exception as e:
            print(f"  Warning: DFT failed for {seq_id}: {e}")
            failed.append(record)

    SeqIO.write(passed, output_fasta, "fasta")

    if scores_file:
        with open(scores_file, "w") as f:
            f.write("sequence_id,homo_eV,lumo_eV,gap_eV\n")
            for seq_id, s in scores.items():
                f.write(f"{seq_id},{s['homo']:.4f},{s['lumo']:.4f},{s['gap']:.4f}\n")
        print(f"  DFT scores saved to {scores_file}")

    total = len(sequences)
    print(f"DFT filter (HOMO-LUMO gap >= {min_homo_lumo_gap} eV):")
    print(f"  Input:  {total} sequences")
    print(f"  Passed: {len(passed)} ({100*len(passed)/total:.1f}%)")
    print(f"  Failed: {len(failed)} sequences")
    print(f"  Output: {output_fasta}")
    return passed


def main():
    parser = argparse.ArgumentParser(description="Filter sequences by DFT electronic properties.")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    parser.add_argument("--structures_dir", required=True, help="Directory with ESMFold PDB structures")
    parser.add_argument("--docking_dir", required=True, help="Directory with docking pose PDB files")
    parser.add_argument("--min_gap", type=float, default=2.0, help="Min HOMO-LUMO gap in eV (default: 2.0)")
    parser.add_argument("--scores_file", default=None, help="Optional: save DFT scores to CSV")
    parser.add_argument("--work_dir", default="dft_work", help="Working directory for DFT files")
    args = parser.parse_args()

    filter_by_dft(
        args.input, args.output, args.structures_dir, args.docking_dir,
        args.min_gap, args.scores_file, args.work_dir
    )


if __name__ == "__main__":
    main()
