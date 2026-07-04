#!/usr/bin/env python3
import sys
import os

def load_lines(path):
    with open(path, 'r') as f:
        return f.readlines()

def write_merged(receptor_lines, ligand_lines, resname, chain, resseq, out_path):
    # Determine last serial used in receptor to continue numbering
    serial = 0
    for line in receptor_lines:
        if line.startswith(('ATOM', 'HETATM')):
            try:
                serial = max(serial, int(line[6:11]))
            except Exception:
                pass
    out = []
    # Write receptor as-is
    for line in receptor_lines:
        if line.startswith('END') or line.startswith('ENDMDL'):
            continue
        out.append(line)
    # Append ligand, rewriting fields
    for line in ligand_lines:
        if not line.startswith(('ATOM', 'HETATM')):
            continue
        name = line[12:16]
        x = line[30:38]
        y = line[38:46]
        z = line[46:54]
        occ = line[54:60] if len(line) >= 60 else '  1.00'
        b = line[60:66] if len(line) >= 66 else '  0.00'
        elem = (line[76:78] if len(line) >= 78 else '  ').rjust(2)
        serial += 1
        new = f"HETATM{serial:5d} {name} {resname:>3} {chain}{resseq:4d}    {x}{y}{z}{occ}{b}          {elem}\n"
        out.append(new)
    out.append('END\n')
    with open(out_path, 'w') as f:
        f.writelines(out)

def parse_mol2(path):
    atoms = []
    bonds = []
    with open(path, 'r') as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.upper().startswith('@<TRIPOS>ATOM'):
            i += 1
            while i < len(lines) and not lines[i].startswith('@'):
                parts = lines[i].split()
                if len(parts) >= 6:
                    name = parts[1][:4]
                    x = float(parts[2])
                    y = float(parts[3])
                    z = float(parts[4])
                    atype = parts[5]
                    elem = ''.join([c for c in atype if c.isalpha()])[:2].strip().title()
                    if not elem:
                        elem = name[0].upper()
                    atoms.append((name, x, y, z, elem))
                i += 1
            continue
        if line.upper().startswith('@<TRIPOS>BOND'):
            i += 1
            while i < len(lines) and not lines[i].startswith('@'):
                parts = lines[i].split()
                if len(parts) >= 3:
                    try:
                        a = int(parts[1])
                        b = int(parts[2])
                        bonds.append((a, b))
                    except Exception:
                        pass
                i += 1
            continue
        i += 1
    return atoms, bonds

def write_merged_with_mol2(receptor_lines, mol2_path, resname, chain, resseq, out_path):
    serial = 0
    for line in receptor_lines:
        if line.startswith(('ATOM', 'HETATM')):
            try:
                serial = max(serial, int(line[6:11]))
            except Exception:
                pass
    out = []
    for line in receptor_lines:
        if line.startswith('END') or line.startswith('ENDMDL'):
            continue
        out.append(line)
    atoms, bonds = parse_mol2(mol2_path)
    idx_to_serial = {}
    for idx, (name, x, y, z, elem) in enumerate(atoms, start=1):
        serial += 1
        idx_to_serial[idx] = serial
        new = (
            f"HETATM{serial:5d} {name:<4} {resname:>3} {chain}{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2}\n"
        )
        out.append(new)
    for a, b in bonds:
        sa = idx_to_serial.get(a)
        sb = idx_to_serial.get(b)
        if sa and sb:
            out.append(f"CONECT{sa:5d}{sb:5d}\n")
    out.append('END\n')
    with open(out_path, 'w') as f:
        f.writelines(out)

if __name__ == '__main__':
    if len(sys.argv) != 7:
        print('Usage: merge_receptor_ligand.py <receptor.pdb> <ligand.pdb> <RESNAME> <CHAIN> <RESSEQ> <out.pdb>')
        sys.exit(1)
    rec, lig, resname, chain, resseq, outp = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]), sys.argv[6]
    if lig.lower().endswith('.mol2'):
        write_merged_with_mol2(load_lines(rec), lig, resname[:3].upper(), chain[:1].upper(), resseq, out_path=outp)
    else:
        write_merged(load_lines(rec), load_lines(lig), resname[:3].upper(), chain[:1].upper(), resseq, out_path=outp)
    
