#!/usr/bin/env python3
import sys
import os
import math
import re

CUTOFF = 4.0  # Angstroms

class Atom:
    def __init__(self, serial, name, resname, chain, resi, x, y, z, element):
        self.serial = serial
        self.name = name.strip()
        self.resname = resname.strip()
        self.chain = chain.strip()
        self.resi = int(resi)
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.element = element.strip()
    def coord(self):
        return (self.x, self.y, self.z)

def parse_pdb(path):
    atoms = []
    with open(path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            try:
                serial = int(line[6:11])
                name = line[12:16]
                resname = line[17:21]
                chain = line[21:22]
                resi = line[22:26]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = (line[76:78].strip() or name.strip()[0])
            except Exception:
                # Fallback: whitespace parsing for atypical alignment
                parts = re.split(r"\s+", line.strip())
                if len(parts) < 9:
                    continue
                try:
                    # HETATM, serial, name, resname, chainResi, x, y, z, ...
                    serial = int(parts[1])
                    name = parts[2].ljust(4)[:4]
                    resname = parts[3].ljust(4)[:4]
                    chainResi = parts[4]
                    m = re.match(r"([A-Za-z]?)(-?\d+)", chainResi)
                    chain = (m.group(1) if m else " ")
                    resi = (m.group(2) if m else "0")
                    x = float(parts[5]); y = float(parts[6]); z = float(parts[7])
                    element = (parts[8] if len(parts) > 8 else name.strip()[0])
                except Exception:
                    continue
            atoms.append(Atom(serial, name, resname, chain, resi, x, y, z, element))
    return atoms

def d(a: Atom, b: Atom):
    dx = a.x - b.x
    dy = a.y - b.y
    dz = a.z - b.z
    return math.sqrt(dx*dx + dy*dy + dz*dz)


def contacts(lig_atoms, other_atoms, cutoff=CUTOFF):
    pairs = []
    mind = float('inf')
    for la in lig_atoms:
        for oa in other_atoms:
            dist = d(la, oa)
            if dist < mind:
                mind = dist
            if dist <= cutoff:
                pairs.append((dist, la, oa))
    pairs.sort(key=lambda x: x[0])
    return mind, pairs


def main():
    if len(sys.argv) < 3:
        print('Usage: compute_contacts.py <merged_complex.pdb> <selector> [out.tsv]')
        sys.exit(1)
    pdb_path = sys.argv[1]
    selector = sys.argv[2]
    out_path = sys.argv[3] if len(sys.argv) > 3 else ''

    atoms = parse_pdb(pdb_path)
    lig = []
    sel = selector.strip()
    if sel.lower().startswith('resn:'):
        lig_resn = sel.split(':',1)[1].strip().upper()
        lig = [a for a in atoms if a.resname.strip().upper() == lig_resn]
    else:
        # chain:Z,resi:3001 form or just resn code fallback
        m_chain = re.search(r"chain:([A-Za-z])", sel, flags=re.I)
        m_resi = re.search(r"resi:(-?\d+)", sel, flags=re.I)
        if m_chain and m_resi:
            ch = m_chain.group(1).upper()
            ri = int(m_resi.group(1))
            lig = [a for a in atoms if a.chain.upper()==ch and a.resi==ri]
        else:
            lig_resn = sel.strip().upper()
            lig = [a for a in atoms if a.resname.strip().upper() == lig_resn]
    fad = [a for a in atoms if a.resname.strip().upper() == 'FAD']
    prot = [a for a in atoms if (a not in lig) and (a.resname.strip().upper() != 'FAD')]

    if not lig:
        print(f'No ligand atoms found for selector {selector} in {pdb_path}')
        sys.exit(2)
    if not fad:
        print('Warning: No FAD found in complex; FAD contacts will be empty')

    min_prot, prot_pairs = contacts(lig, prot)
    min_fad, fad_pairs = contacts(lig, fad) if fad else (float('inf'), [])

    lines = []
    lines.append(f'# Complex: {os.path.basename(pdb_path)}\n')
    lines.append(f'# Selector: {selector}\n')
    lines.append(f'MIN_DIST_LIG_PROT\t{min_prot:.3f}\n')
    if fad:
        lines.append(f'MIN_DIST_LIG_FAD\t{min_fad:.3f}\n')
    lines.append(f'NUM_CONTACTS_LIG_PROT_{CUTOFF:.1f}A\t{len(prot_pairs)}\n')
    if fad:
        lines.append(f'NUM_CONTACTS_LIG_FAD_{CUTOFF:.1f}A\t{len(fad_pairs)}\n')
    lines.append('# Top 15 ligand-protein contacts (dist, lig_atom, prot_atom, prot_res)\n')
    for dist, la, oa in prot_pairs[:15]:
        lines.append(f'{dist:.3f}\t{la.name}@{la.resname}{la.resi}\t{oa.name}@{oa.resname}{oa.resi}\t{oa.chain}\n')
    if fad_pairs:
        lines.append('# Top 15 ligand-FAD contacts (dist, lig_atom, fad_atom)\n')
        for dist, la, oa in fad_pairs[:15]:
            lines.append(f'{dist:.3f}\t{la.name}@{la.resname}{la.resi}\t{oa.name}@{oa.resname}{oa.resi}\n')

    if out_path:
        with open(out_path, 'w') as f:
            f.writelines(lines)
        print(f'Wrote {out_path}')
    else:
        sys.stdout.writelines(lines)

if __name__ == '__main__':
    main()
