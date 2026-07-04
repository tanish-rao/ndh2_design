#!/usr/bin/env python3
import argparse
import os
import re
from pathlib import Path
from typing import Dict, List
from Bio import SeqIO

IDX_RE = re.compile(r".*?(\d+)$")

def parse_index(seq_id: str) -> int:
    m = IDX_RE.match(seq_id.replace('>', '').strip())
    if not m:
        raise ValueError(f"Cannot parse index from id: {seq_id}")
    return int(m.group(1))


def load_codons(source_codons: Path) -> List[str]:
    with open(source_codons) as f:
        return [ln.strip() for ln in f if ln.strip()]


def write_subset_for_fasta(fasta_path: Path, codon_lines: List[str], out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    # Copy proteins to tier folder
    prot_out = out_dir / 'proteins.fasta'
    # To maintain order, re-parse and write
    records = list(SeqIO.parse(str(fasta_path), 'fasta'))
    SeqIO.write(records, str(prot_out), 'fasta')
    # Write matching codon lines, based on numeric index in header
    codon_out = out_dir / 'codons.txt'
    with open(codon_out, 'w') as fout:
        for rec in records:
            idx = parse_index(rec.id)
            if idx < 0 or idx >= len(codon_lines):
                continue
            fout.write(codon_lines[idx] + '\n')
    return prot_out, codon_out


def discover_tier_fastas(results_dir: Path) -> Dict[int, Path]:
    mapping = {}
    for tier in range(1, 8):
        # canonical name written by orchestrator
        canon = results_dir / {
            1: 'tier1_length.fasta',
            2: 'tier2_blast.fasta',
            3: 'tier3_esmfold.fasta',
            4: 'tier4_topology.fasta',
            5: 'tier5_clustered.fasta',
            6: 'tier6_docking.fasta',
            7: 'tier7_dft.fasta',
        }[tier]
        if canon.exists():
            mapping[tier] = canon
    return mapping


def main():
    ap = argparse.ArgumentParser(description='Create per-tier codon/protein outputs')
    ap.add_argument('--results_dir', required=True, help='Pipeline results dir (contains tier*.fasta)')
    ap.add_argument('--source_codons', required=True, help='Codon tokens used to create Tier 0 proteins')
    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    codon_lines = load_codons(Path(args.source_codons))
    tiers = discover_tier_fastas(results_dir)

    for tier, fasta in sorted(tiers.items()):
        tier_dir = results_dir / f'tier{tier}'
        prot_out, codon_out = write_subset_for_fasta(fasta, codon_lines, tier_dir)
        print(f'Tier {tier}: wrote {prot_out} and {codon_out}')

if __name__ == '__main__':
    main()
