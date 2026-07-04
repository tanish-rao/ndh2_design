#!/usr/bin/env python3
import os
import numpy as np

WORKDIR = os.path.abspath('outputs/fad_refine_64538480')
PRE = os.path.join(WORKDIR, 'fad_out.pdb')
POST = os.path.join(WORKDIR, 'fad_minimized.pdb')


def read_pdb_coords(path):
    atoms = []
    with open(path, 'r') as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            # Skip hydrogens
            element = (line[76:78].strip() or line[12:16].strip()[0]).upper()
            if element == 'H':
                continue
            name = line[12:16].strip()
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            atoms.append((name, np.array([x, y, z], dtype=float)))
    return atoms


def kabsch(P, Q):
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
    pre_atoms = read_pdb_coords(PRE)
    post_atoms = read_pdb_coords(POST)

    # Match by atom name order intersection
    pre_map = {a[0]: a[1] for a in pre_atoms}
    post_map = {a[0]: a[1] for a in post_atoms}
    names = [n for n in pre_map.keys() if n in post_map]
    if len(names) < 5:
        raise SystemExit('Not enough matching heavy atoms between pre and post to compute RMSD.')

    P = np.stack([pre_map[n] for n in names], axis=0)
    Q = np.stack([post_map[n] for n in names], axis=0)

    # Align pre to post
    R, t = kabsch(P, Q)
    P_aln = (R @ P.T).T + t

    diffs = P_aln - Q
    rmsd = np.sqrt((diffs**2).sum(axis=1).mean())

    print(f"Matched atoms: {len(names)}")
    print(f"Heavy-atom RMSD (Å): {rmsd:.3f}")

if __name__ == '__main__':
    main()
