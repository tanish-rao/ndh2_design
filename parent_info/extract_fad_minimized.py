#!/usr/bin/env python3
import os

WORKDIR = os.path.abspath('outputs/fad_refine_64538480')
INP = os.path.join(WORKDIR, 'complex_receptor_fad_minimized.pdb')
OUT = os.path.join(WORKDIR, 'fad_minimized.pdb')

resnames = {"FAD", "FAD0", "FMN", "FADH", "FADH2"}

with open(INP, 'r') as fin, open(OUT, 'w') as fout:
    serial = 1
    for line in fin:
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue
        resname = line[17:20].strip()
        if resname in resnames or resname.upper().startswith('FAD'):
            # Heavy atoms only
            element = (line[76:78].strip() or line[12:16].strip()[0]).upper()
            if element == 'H':
                continue
            # Re-emit with continuous serials
            new = f"{line[:6]}{serial:5d}{line[11:]}"
            fout.write(new)
            serial += 1
    fout.write('END\n')
print('Wrote', OUT)
