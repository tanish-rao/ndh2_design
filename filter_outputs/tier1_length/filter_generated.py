#!/usr/bin/env python3
import os
import sys
import argparse
import glob
import json
from datetime import datetime
from typing import List, Dict, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

STOP_CODONS = {"TAA", "TAG", "TGA"}


def gc_content(codons: List[str]) -> float:
    seq = ''.join(codons).upper()
    if not seq:
        return 0.0
    g = seq.count('G')
    c = seq.count('C')
    a = seq.count('A')
    t = seq.count('T')
    total = g + c + a + t
    return ((g + c) / total * 100.0) if total else 0.0


def trim_trailing_stops(codons: List[str]) -> Tuple[List[str], int]:
    n = 0
    while codons and codons[-1].upper() in STOP_CODONS:
        codons = codons[:-1]
        n += 1
    return codons, n


def has_internal_stop(codons: List[str]) -> bool:
    if not codons:
        return False
    for c in codons[:-1]:
        if c.upper() in STOP_CODONS:
            return True
    return False


def load_inputs(files: List[str], pattern: str) -> List[str]:
    paths = []
    for f in files:
        if f:
            paths.append(f)
    if pattern:
        paths.extend(sorted(glob.glob(pattern)))
    # Default: all generated files in default directory
    if not paths:
        paths = sorted(glob.glob('train_model/generation/generated_codons_*.txt'))
    # Deduplicate while preserving order
    seen = set()
    uniq = []
    for p in paths:
        if p not in seen and os.path.isfile(p):
            seen.add(p)
            uniq.append(p)
    return uniq


def filter_and_summarize(inputs: List[str], out_dir: str, min_aa: int = 330, max_aa: int = 600) -> Dict:
    os.makedirs(out_dir, exist_ok=True)

    kept: List[str] = []
    kept_meta: List[Dict] = []
    dropped: List[Dict] = []

    seen_trimmed = set()

    total_in = 0

    for path in inputs:
        with open(path, 'r') as f:
            lines = [ln.strip() for ln in f if ln.strip()]
        for i, line in enumerate(lines):
            total_in += 1
            toks = line.split()
            toks_trim, n_trim = trim_trailing_stops(toks)
            aa_len = len(toks_trim)
            internal = has_internal_stop(toks)
            reason = None
            if internal:
                reason = 'internal_stop'
            elif not (min_aa <= aa_len <= max_aa):
                reason = 'length_out_of_range'
            elif ' '.join(toks_trim) in seen_trimmed:
                reason = 'duplicate_after_trim'

            if reason is None:
                seen_trimmed.add(' '.join(toks_trim))
                gc = gc_content(toks_trim)
                kept.append(' '.join(toks_trim))
                kept_meta.append({
                    'source_file': os.path.abspath(path),
                    'source_index': i,
                    'aa_len': aa_len,
                    'gc_pct': gc,
                    'trimmed_trailing_stops': n_trim,
                })
            else:
                dropped.append({
                    'source_file': os.path.abspath(path),
                    'source_index': i,
                    'reason': reason,
                    'orig_len_codons': len(toks),
                    'trimmed_len_codons': aa_len,
                })

    # Write outputs
    cleaned_path = os.path.join(out_dir, 'generated_cleaned_codons.txt')
    with open(cleaned_path, 'w') as f:
        for s in kept:
            f.write(s + '\n')

    with open(os.path.join(out_dir, 'dropped.tsv'), 'w') as f:
        if dropped:
            headers = list(dropped[0].keys())
            f.write('\t'.join(headers) + '\n')
            for r in dropped:
                f.write('\t'.join(str(r[h]) for h in headers) + '\n')
        else:
            f.write('source_file\tsource_index\treason\torig_len_codons\ttrimmed_len_codons\n')

    meta_path = os.path.join(out_dir, 'kept_metadata.tsv')
    with open(meta_path, 'w') as f:
        if kept_meta:
            headers = list(kept_meta[0].keys())
            f.write('\t'.join(headers) + '\n')
            for r in kept_meta:
                f.write('\t'.join(str(r[h]) for h in headers) + '\n')
        else:
            f.write('source_file\tsource_index\taa_len\tgc_pct\ttrimmed_trailing_stops\n')

    # Summary & plots
    lengths = [r['aa_len'] for r in kept_meta]
    gcs = [r['gc_pct'] for r in kept_meta]

    summary = {
        'inputs': [os.path.abspath(p) for p in inputs],
        'output_cleaned': os.path.abspath(cleaned_path),
        'n_input_total': total_in,
        'n_kept': len(kept),
        'n_dropped': len(dropped),
        'pct_kept': (len(kept) / total_in * 100.0) if total_in else 0.0,
        'min_aa': min_aa,
        'max_aa': max_aa,
    }

    with open(os.path.join(out_dir, 'summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)

    plt.figure(figsize=(6,4))
    plt.hist(lengths, bins=40, color='#4e79a7')
    plt.xlabel('AA length (after trimming stops)')
    plt.ylabel('Count')
    plt.title('Kept AA length distribution')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'kept_hist_aa_length.png'), dpi=150)
    plt.close()

    plt.figure(figsize=(6,4))
    plt.hist(gcs, bins=40, color='#59a14f')
    plt.xlabel('GC%')
    plt.ylabel('Count')
    plt.title('Kept GC% distribution')
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'kept_hist_gc_pct.png'), dpi=150)
    plt.close()

    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inputs', nargs='*', default=[], help='List of input codon files (one per line)')
    ap.add_argument('--glob', dest='pattern', default='', help='Glob pattern for inputs (e.g., train_model/generation/generated_codons_*.txt)')
    ap.add_argument('--out_dir', default=None, help='Output directory (default: train_model/generation/filtered_<timestamp>)')
    ap.add_argument('--min_aa', type=int, default=330)
    ap.add_argument('--max_aa', type=int, default=600)
    args = ap.parse_args()

    inputs = load_inputs(args.inputs, args.pattern)
    if not inputs:
        print('No input files found. Specify --inputs or --glob pattern.', file=sys.stderr)
        sys.exit(1)

    if args.out_dir is None:
        ts = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.out_dir = os.path.join('train_model', 'generation', f'filtered_{ts}')

    summary = filter_and_summarize(inputs, args.out_dir, args.min_aa, args.max_aa)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
