#!/usr/bin/env python3
"""
Energy-minimize the NDH2 + FAD complex using OpenMM.
- Protein: Amber14 protein force field (implicit solvent OBC2)
- Ligand (FAD): SMIRNOFF via OpenFF (requires RDKit + openff-toolkit + openmmforcefields)

Inputs (defaults target the recent redock output directory):
  --workdir PATH            Directory containing receptor.pdb and fad_out.pdb (default: outputs/fad_refine_64538480)
  --protein PDB             Protein PDB path (default: {workdir}/receptor.pdb)
  --ligand PDB              Ligand PDB path (default: {workdir}/fad_out.pdb)
  --outdir PATH             Output directory (default: {workdir})

Outputs:
  - complex_receptor_fad_minimized.pdb
  - openmm_minimize_energy.txt (final potential energy, kJ/mol)
"""
import argparse
import os
import subprocess
import shutil

from openmm import unit
from openmm.app import PDBFile, Modeller, Simulation
from openmm.app import ForceField
from openmm.app import NoCutoff, HBonds
from openmm import LangevinIntegrator, Platform, LocalEnergyMinimizer, Vec3
from openmm.app import AmberPrmtopFile, AmberInpcrdFile

from openff.toolkit.topology import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from rdkit import Chem
from rdkit.Chem import AllChem


def _ensure_sdf_from_ligand(ligand_pdb: str) -> str:
    base, _ = os.path.splitext(ligand_pdb)
    sdf = base + '.sdf'
    if os.path.exists(sdf):
        return sdf
    # Try to convert from a PDBQT with same basename if present (better bond orders)
    pdbqt = base + '.pdbqt'
    obabel = shutil.which('obabel') or shutil.which('obabel.exe')
    if obabel:
        if os.path.exists(pdbqt):
            subprocess.run([obabel, '-ipdbqt', pdbqt, '-osdf', '-O', sdf], check=True)
            if os.path.exists(sdf):
                return sdf
        # fallback: convert from PDB
        subprocess.run([obabel, '-ipdb', ligand_pdb, '-osdf', '-O', sdf], check=True)
        if os.path.exists(sdf):
            return sdf
    return ''


def build_system(protein_pdb: str, ligand_pdb: str):
    # Load protein and ligand
    protein = PDBFile(protein_pdb)
    # Prefer chemistry from SDF (with proper bond orders). Auto-generate via Open Babel if available.
    sdf = _ensure_sdf_from_ligand(ligand_pdb)
    lig_mol = None
    if sdf:
        lig_mol = Molecule.from_file(sdf)
    else:
        # Fallback: build RDKit molecule from PDB and preserve 3D conformer
        rdk = Chem.MolFromPDBFile(ligand_pdb, removeHs=False)
        if rdk is None:
            raise RuntimeError(f"Failed to determine ligand chemistry from {ligand_pdb}. Consider installing Open Babel to enable SDF generation.")
        if rdk.GetNumConformers() == 0:
            AllChem.EmbedMolecule(rdk, randomSeed=0xF00D)
        lig_mol = Molecule.from_rdkit(rdk, allow_undefined_stereo=True)
    # Build OpenMM topology/positions for ligand from OpenFF molecule conformer
    if not lig_mol.conformers:
        raise RuntimeError("Ligand molecule lacks 3D coordinates; ensure SDF generation succeeded.")
    lig_top_omm = lig_mol.to_topology().to_openmm()
    # Convert conformer to list of Vec3 with nanometer units
    if hasattr(lig_mol.conformers[0], 'value_in_unit'):
        conf_nm = lig_mol.conformers[0].value_in_unit(unit.nanometer)
    else:
        # Assume Angstroms if unitless numpy; convert to nm
        import numpy as np
        arr = np.array(lig_mol.conformers[0], dtype=float) / 10.0
        conf_nm = arr
    lig_pos_omm = [Vec3(float(x), float(y), float(z)) * unit.nanometer for x, y, z in conf_nm]

    modeller = Modeller(protein.topology, protein.positions)
    modeller.add(lig_top_omm, lig_pos_omm)

    # Add hydrogens at pH 7.0 (protein and ligand)
    modeller.addHydrogens(pH=7.0)

    # Prepare force field: Amber14 protein + implicit solvent; SMIRNOFF for ligand
    smirnoff = SMIRNOFFTemplateGenerator(molecules=[lig_mol],
                                         forcefield='openff-2.0.0.offxml')
    ff = ForceField('amber14/protein.ff14SB.xml', 'implicit/obc2.xml')
    ff.registerTemplateGenerator(smirnoff.generator)

    # Create system; constraints on bonds to H
    system = ff.createSystem(modeller.topology,
                             nonbondedMethod=NoCutoff,
                             constraints=HBonds)
    return system, modeller


def minimize(protein_pdb: str, ligand_pdb: str, out_pdb: str, energy_txt: str):
    # If Amber parameter files are present in the ligand directory, use them directly.
    workdir = os.path.dirname(ligand_pdb)
    prmtop_path = os.path.join(workdir, 'complex.prmtop')
    inpcrd_path = os.path.join(workdir, 'complex.inpcrd')

    if os.path.exists(prmtop_path) and os.path.exists(inpcrd_path):
        prmtop = AmberPrmtopFile(prmtop_path)
        inpcrd = AmberInpcrdFile(inpcrd_path)
        system = prmtop.createSystem(nonbondedMethod=NoCutoff, constraints=HBonds)
        integrator = LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.004*unit.picoseconds)
        platform = Platform.getPlatformByName('CPU')
        sim = Simulation(prmtop.topology, system, integrator, platform)
        sim.context.setPositions(inpcrd.positions)
    else:
        system, modeller = build_system(protein_pdb, ligand_pdb)
        integrator = LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.004*unit.picoseconds)
        platform = Platform.getPlatformByName('CPU')
        sim = Simulation(modeller.topology, system, integrator, platform)
        sim.context.setPositions(modeller.positions)

    # Energy before
    state0 = sim.context.getState(getEnergy=True)
    e0 = state0.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    # Minimize
    LocalEnergyMinimizer.minimize(sim.context, tolerance=10.0*unit.kilojoule_per_mole/unit.nanometer, maxIterations=5000)

    # Energy after
    state = sim.context.getState(getPositions=True, getEnergy=True)
    e1 = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    # Save minimized structure
    with open(out_pdb, 'w') as f:
        PDBFile.writeFile(sim.topology, state.getPositions(), f, keepIds=True)

    with open(energy_txt, 'w') as f:
        f.write(f"Initial potential energy (kJ/mol): {e0:.2f}\n")
        f.write(f"Final potential energy   (kJ/mol): {e1:.2f}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workdir', default='outputs/fad_refine_64538480')
    ap.add_argument('--protein', default=None)
    ap.add_argument('--ligand', default=None)
    ap.add_argument('--outdir', default=None)
    args = ap.parse_args()

    workdir = os.path.abspath(args.workdir)
    protein = args.protein or os.path.join(workdir, 'receptor.pdb')
    ligand = args.ligand or os.path.join(workdir, 'fad_out.pdb')
    outdir = os.path.abspath(args.outdir or workdir)
    os.makedirs(outdir, exist_ok=True)

    out_pdb = os.path.join(outdir, 'complex_receptor_fad_minimized.pdb')
    energy_txt = os.path.join(outdir, 'openmm_minimize_energy.txt')

    print('Protein PDB:', protein)
    print('Ligand  PDB:', ligand)
    print('Out dir    :', outdir)

    minimize(protein, ligand, out_pdb, energy_txt)
    print('Minimization complete.')
    print('Wrote:', out_pdb)
    print('Energy report:', energy_txt)


if __name__ == '__main__':
    main()
