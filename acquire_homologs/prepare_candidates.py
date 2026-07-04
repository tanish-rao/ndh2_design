#!/usr/bin/env python3
import os
import sys
import re
from datetime import datetime

def ts():
    return datetime.now().isoformat(timespec='seconds')

GENETIC_CODE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

EXCLUDE_TERMS = [
    'nuod','nuoc','nuob','nuoh','nuoi','nuoj','nuok','nuol','nuom','nuon',
    'complex i','nadh-ubiquinone oxidoreductase','nadh:ubiquinone',
    'nadh dehydrogenase subunit','nadh dehydrogenase i',
]

def is_excluded(header: str) -> bool:
    hl = header.lower()
    for t in EXCLUDE_TERMS:
        if t in hl:
            return True
    return False

def translate(seq: str) -> str:
    aa = []
    for i in range(0, len(seq)-2, 3):
        aa.append(GENETIC_CODE.get(seq[i:i+3], 'X'))
    return ''.join(aa)

def fasta_stream(path):
    with open(path) as f:
        header = None
        parts = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(parts).upper().replace(' ', '')
                header = line[1:].strip()
                parts = []
            else:
                parts.append(line)
        if header is not None:
            yield header, ''.join(parts).upper().replace(' ', '')

def main():
    # Args: batch_dir out_dir [min_aa] [max_aa]
    if len(sys.argv) < 3:
        print('usage: prepare_candidates.py <batch_dir> <out_dir> [min_aa] [max_aa]')
        sys.exit(1)
    batch_dir = sys.argv[1]
    out_dir = sys.argv[2]
    min_aa = int(sys.argv[3]) if len(sys.argv) > 3 else 330
    max_aa = int(sys.argv[4]) if len(sys.argv) > 4 else 600
    os.makedirs(out_dir, exist_ok=True)
    cand_faa = os.path.join(out_dir, 'candidates.faa')
    map_tsv = os.path.join(out_dir, 'candidates_map.tsv')

    batch_files = sorted([os.path.join(batch_dir, x) for x in os.listdir(batch_dir) if re.match(r'^batch_\d+\.fa$', x)])
    total = 0
    kept = 0
    with open(cand_faa, 'w') as fout, open(map_tsv, 'w') as fmap:
        fmap.write('name\tprotein_id\tbatch_file\torig_header\n')
        for p in batch_files:
            for h, s in fasta_stream(p):
                total += 1
                if is_excluded(h):
                    continue
                aa = translate(s)
                if len(aa) < min_aa or len(aa) > max_aa:
                    continue
                if aa.count('X') > 5:
                    continue
                # get protein_id if present
                m = re.search(r'protein_id=([^\] ]+)', h)
                pid = m.group(1) if m else None
                name = pid if pid else f'noid_{abs(hash(h)) % (10**10)}'
                fout.write(f'>{name}\n')
                for i in range(0, len(aa), 60):
                    fout.write(aa[i:i+60] + '\n')
                fmap.write(f'{name}\t{pid or ""}\t{os.path.basename(p)}\t{h}\n')
                kept += 1
    print(f'[{ts()}] Wrote candidates: {kept}/{total} kept -> {cand_faa}')

if __name__ == '__main__':
    main()
