#!/usr/bin/env python3
import os
import sys
import json
import argparse
from datetime import datetime
from typing import List, Tuple, Dict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

STOP_CODONS = {"TAA", "TAG", "TGA"}


def gc_content(codons: List[str]) -> Tuple[int, int, float]:
    seq = ''.join(codons)
    seq = seq.upper()
    g = seq.count('G')
    c = seq.count('C')
    a = seq.count('A')
    t = seq.count('T')
    total = g + c + a + t
    gc_pct = (g + c) / total * 100.0 if total > 0 else 0.0
    return g + c, a + t, gc_pct


def trim_trailing_stops(codons: List[str]) -> Tuple[List[str], int]:
    n = 0
    while codons and codons[-1].upper() in STOP_CODONS:
        codons = codons[:-1]
        n += 1
    return codons, n


def has_internal_stop(codons: List[str]) -> bool:
    if not codons:
        return False
    # Allow trailing stops (handled separately). Check all but last.
    for i, c in enumerate(codons[:-1]):
        if c.upper() in STOP_CODONS:
            return True
    return False


def analyze_file(codon_file: str, out_dir: str):
    os.makedirs(out_dir, exist_ok=True)

    per_seq: List[Dict] = []
    unique_set = set()
    duplicate_idx = []

    with open(codon_file, 'r') as f:
        lines = [ln.strip() for ln in f if ln.strip()]

    for i, line in enumerate(lines):
        toks = line.split()
        toks_trimmed, n_trim = trim_trailing_stops(toks)
        aa_len = len(toks_trimmed)
        gc, at, gc_pct = gc_content(toks_trimmed)
        internal_stop = has_internal_stop(toks)
        key = ' '.join(toks_trimmed)
        is_dup = key in unique_set
        if not is_dup:
            unique_set.add(key)
        else:
            duplicate_idx.append(i)
        per_seq.append({
            'index': i,
            'n_codons': len(toks),
            'aa_len': aa_len,
            'trimmed_trailing_stops': n_trim,
            'gc_pct': gc_pct,
            'gc_bases': gc,
            'at_bases': at,
            'has_internal_stop': internal_stop,
            'is_duplicate': is_dup,
            'in_len_range_330_600': 330 <= aa_len <= 600,
        })

    # Aggregate stats
    total = len(per_seq)
    unique = len(unique_set)
    in_range = sum(1 for r in per_seq if r['in_len_range_330_600'])
    dups = total - unique
    internal_stop_n = sum(1 for r in per_seq if r['has_internal_stop'])

    summary = {
        'input_file': os.path.abspath(codon_file),
        'total_sequences': total,
        'unique_after_trimming': unique,
        'duplicates_after_trimming': dups,
        'in_len_range_330_600': in_range,
        'pct_in_len_range': (in_range / total * 100.0) if total else 0.0,
        'sequences_with_internal_stop': internal_stop_n,
    }

    # Write TSV and JSON
    tsv_path = os.path.join(out_dir, 'generated_metrics.tsv')
    with open(tsv_path, 'w') as f:
        header = ['index','n_codons','aa_len','trimmed_trailing_stops','gc_pct','gc_bases','at_bases','has_internal_stop','is_duplicate','in_len_range_330_600']
        f.write('\t'.join(header) + '\n')
        for r in per_seq:
            f.write('\t'.join(str(r[h]) for h in header) + '\n')

    with open(os.path.join(out_dir, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    # Plots
    lengths = [r['aa_len'] for r in per_seq]
    gcs = [r['gc_pct'] for r in per_seq]

    plt.figure(figsize=(6,4))
    plt.hist(lengths, bins=40, color='#4e79a7')
    plt.xlabel('AA length (codons after trimming stops)')
    plt.ylabel('Count')
    plt.title('Generated AA length distribution')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'hist_aa_length.png'), dpi=150)
    plt.close()

    plt.figure(figsize=(6,4))
    plt.hist(gcs, bins=40, color='#59a14f')
    plt.xlabel('GC%')
    plt.ylabel('Count')
    plt.title('Generated GC% distribution')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'hist_gc_pct.png'), dpi=150)
    plt.close()

    plt.figure(figsize=(5,5))
    plt.scatter(lengths, gcs, s=6, alpha=0.4, color='#e15759')
    plt.xlabel('AA length')
    plt.ylabel('GC%')
    plt.title('Length vs GC%')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'len_vs_gc.png'), dpi=150)
    plt.close()

    print(json.dumps(summary, indent=2))
    print(f"Wrote: {tsv_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--codon_file', required=True, help='Path to generated codon tokens (one sequence per line)')
    ap.add_argument('--out_dir', default=None, help='Output directory for diagnostics')
    args = ap.parse_args()

    if args.out_dir is None:
        ts = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.out_dir = os.path.join('train_model', 'generation', f'diag_{ts}')

    analyze_file(args.codon_file, args.out_dir)


if __name__ == '__main__':
    main()
