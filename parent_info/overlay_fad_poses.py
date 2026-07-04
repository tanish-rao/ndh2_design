#!/usr/bin/env python3
import os
import math
import numpy as np

POST_DIR = os.path.abspath('outputs/fad_refine_64538480')
PRE_RECEPTOR = os.path.join(POST_DIR, 'receptor.pdb')
PRE_FAD = os.path.join(POST_DIR, 'fad_out.pdb')
POST_COMPLEX = os.path.join(POST_DIR, 'complex_receptor_fad_minimized.pdb')
OUT_PDB = os.path.join(POST_DIR, 'overlay_fad_pre_vs_post.pdb')


def parse_pdb_atoms(path):
    atoms = []
    with open(path, 'r') as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            rec = line[0:6].strip()
            serial = int(line[6:11]) if line[6:11].strip() else None
            name = line[12:16].strip()
            alt = line[16].strip()
            resname = line[17:20].strip()
            chain = line[21].strip() or 'A'
            resseq = int(line[22:26]) if line[22:26].strip() else 0
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            occ = line[54:60].strip() or '1.00'
            b = line[60:66].strip() or '0.00'
            element = (line[76:78].strip() or name.strip()[0]).upper()
            atoms.append({
                'rec': rec,
                'serial': serial,
                'name': name,
                'alt': alt,
                'resname': resname,
                'chain': chain,
                'resseq': resseq,
                'x': x, 'y': y, 'z': z,
                'occ': occ,
                'b': b,
                'element': element,
                'raw': line.rstrip('\n'),
            })
    return atoms


def write_pdb_atoms(atoms, other_lines, out_path):
    serial = 1
    with open(out_path, 'w') as f:
        for line in other_lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Re-emit with continuous serials
                rec = line[0:6]
                new = f"{rec}{serial:5d}" + line[11:]
                f.write(new)
                serial += 1
            else:
                f.write(line)
        # Append transformed pre-min FAD as distinct residue FAD0 on chain Z
        for a in atoms:
            rec = 'HETATM'
            name = a['name'][:4].rjust(4)
            resname = 'FAD0'
            chain = 'Z'
            resseq = 2001
            x, y, z = a['x'], a['y'], a['z']
            occ = '1.00'
            b = '0.00'
            element = (a['element'] or a['name'][0]).rjust(2)
            line = (
                f"{rec:<6}{serial:5d} {name} {resname:>3} {chain}{resseq:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{float(occ):6.2f}{float(b):6.2f}          {element}\n"
            )
            f.write(line)
            serial += 1
        f.write('END\n')


def kabsch(P, Q):
    # P, Q are Nx3
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    t = Q.mean(axis=0) - R @ P.mean(axis=0)
    return R, t


def main():
    # Load post-minimized complex (reference frame)
    with open(POST_COMPLEX, 'r') as f:
        post_lines = f.readlines()
    post_atoms = parse_pdb_atoms(POST_COMPLEX)
    post_ca = np.array([[a['x'], a['y'], a['z']] for a in post_atoms if a['rec']=='ATOM' and a['name']=='CA'])

    # Load pre-min receptor and FAD
    pre_receptor_atoms = parse_pdb_atoms(PRE_RECEPTOR)
    pre_ca = np.array([[a['x'], a['y'], a['z']] for a in pre_receptor_atoms if a['rec']=='ATOM' and a['name']=='CA'])
    fad_atoms = [a for a in parse_pdb_atoms(PRE_FAD) if a['rec'] in ('ATOM','HETATM')]

    # Align pre-min receptor to post-min complex using CA atoms
    n = min(len(pre_ca), len(post_ca))
    if n < 3:
        raise RuntimeError('Not enough CA atoms to align.')
    R, t = kabsch(pre_ca[:n], post_ca[:n])

    # Transform FAD coordinates into post-min frame
    fad_transformed = []
    for a in fad_atoms:
        v = np.array([a['x'], a['y'], a['z']])
        vt = R @ v + t
        b = a.copy()
        b['x'], b['y'], b['z'] = vt.tolist()
        fad_transformed.append(b)

    # Write overlay: post-min complex + transformed pre-min FAD labeled FAD0 chain Z
    write_pdb_atoms(fad_transformed, post_lines, OUT_PDB)
    print('Wrote overlay:', OUT_PDB)

if __name__ == '__main__':
    main()
