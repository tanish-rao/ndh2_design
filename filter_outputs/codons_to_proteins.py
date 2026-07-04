#!/usr/bin/env python3
import argparse
import os
from typing import List

CODON_TABLE = {
    # Phenylalanine
    'TTT':'F','TTC':'F',
    # Leucine
    'TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    # Isoleucine
    'ATT':'I','ATC':'I','ATA':'I',
    # Methionine (start)
    'ATG':'M',
    # Valine
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    # Serine
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','AGT':'S','AGC':'S',
    # Proline
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    # Threonine
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    # Alanine
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    # Tyrosine
    'TAT':'Y','TAC':'Y',
    # Histidine
    'CAT':'H','CAC':'H',
    # Glutamine
    'CAA':'Q','CAG':'Q',
    # Asparagine
    'AAT':'N','AAC':'N',
    # Lysine
    'AAA':'K','AAG':'K',
    # Aspartic Acid
    'GAT':'D','GAC':'D',
    # Glutamic Acid
    'GAA':'E','GAG':'E',
    # Cysteine
    'TGT':'C','TGC':'C',
    # Tryptophan
    'TGG':'W',
    # Arginine
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGA':'R','AGG':'R',
    # Glycine
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
STOPS = {'TAA','TAG','TGA'}

def translate_codons(codons: List[str]) -> str:
    aas = []
    for c in codons:
        u = c.upper()
        if u in STOPS:
            break
        aa = CODON_TABLE.get(u)
        if aa is None:
            # Unknown codon; skip this sequence entirely
            return ''
        aas.append(aa)
    return ''.join(aas)


def convert_file(input_codons: str, output_fasta: str, prefix: str = 'gen') -> int:
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
    n = 0
    with open(input_codons) as fin, open(output_fasta, 'w') as fout:
        for i, line in enumerate(fin):
            line = line.strip()
            if not line:
                continue
            codons = line.split()
            prot = translate_codons(codons)
            if not prot:
                continue
            n += 1
            header = f'>{prefix}_{i}\n'
            seq = prot + '\n'
            fout.write(header)
            # wrap to 60 columns
            for j in range(0, len(prot), 60):
                fout.write(prot[j:j+60] + '\n')
    return n


def main():
    ap = argparse.ArgumentParser(description='Convert codon-token sequences to protein FASTA')
    ap.add_argument('--input', required=True, help='Path to cleaned codon tokens file')
    ap.add_argument('--output', required=True, help='Output protein FASTA path')
    ap.add_argument('--prefix', default='gen', help='Header prefix')
    args = ap.parse_args()
    n = convert_file(args.input, args.output, args.prefix)
    print(f'Wrote {n} protein sequences to {args.output}')

if __name__ == '__main__':
    main()
