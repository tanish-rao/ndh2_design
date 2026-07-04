#!/usr/bin/env python3
import os
import sys
import re
import statistics as stats
from collections import Counter, defaultdict

FASTA_PATH = 'acquire_homologs/outputs/homologs_filtered.fasta'
META_PATH = 'acquire_homologs/outputs/homologs_metadata.tsv'
OUT_DIR = 'acquire_homologs/outputs'

CODONS = [a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG']

def load_fasta(path):
    entries = []
    header = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header:
                    entries.append((header, ''.join(seq).upper()))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header:
        entries.append((header, ''.join(seq).upper()))
    return entries

def parse_meta(path):
    rows = []
    with open(path) as f:
        next(f)  # header
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            h, nt_len, aa_len = parts[0], int(parts[1]), int(parts[2])
            rows.append((h, nt_len, aa_len))
    return rows

def extract_organism(header: str) -> str:
    # Try to extract [organism=...] tag if present
    m = re.search(r"\[organism=([^\]]+)\]", header)
    if m:
        return m.group(1)
    # Fallback: try to capture between spaces followed by 'gene=' or similar tags
    # or split by hypothetical formats, use first two tokens
    toks = re.split(r"\s+", header)
    if toks:
        return ' '.join(toks[:2])
    return 'unknown'

def gc_content(seq: str) -> float:
    g = seq.count('G')
    c = seq.count('C')
    return 100.0 * (g + c) / max(1, len(seq))

def codon_usage(seq: str) -> Counter:
    cnt = Counter()
    L = (len(seq)//3)*3
    for i in range(0, L, 3):
        cod = seq[i:i+3]
        if len(cod)==3 and all(ch in 'ATGC' for ch in cod):
            cnt[cod] += 1
    return cnt

def write_tsv(path, rows, header=None):
    with open(path, 'w') as f:
        if header:
            f.write('\t'.join(header) + '\n')
        for r in rows:
            f.write('\t'.join(map(str, r)) + '\n')

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    fasta = load_fasta(FASTA_PATH)
    meta = {h: (nt, aa) for h, nt, aa in parse_meta(META_PATH)}

    lengths = []
    gcs = []
    org_counts = Counter()
    agg_codon = Counter()

    per_seq_rows = []

    for header, seq in fasta:
        nt_len, aa_len = meta.get(header, (len(seq), len(seq)//3))
        lengths.append(nt_len)
        gc = gc_content(seq)
        gcs.append(gc)
        org = extract_organism(header)
        org_counts[org] += 1
        cu = codon_usage(seq)
        agg_codon.update(cu)
        per_seq_rows.append((header, org, nt_len, aa_len, f"{gc:.2f}"))

    # Summaries
    summary_lines = []
    summary_lines.append(f"Total sequences: {len(fasta)}\n")
    if lengths:
        summary_lines.append(f"Length nt: mean {stats.mean(lengths):.1f}, median {stats.median(lengths):.1f}, min {min(lengths)}, max {max(lengths)}\n")
    if gcs:
        summary_lines.append(f"GC%: mean {stats.mean(gcs):.2f}, median {stats.median(gcs):.2f}, min {min(gcs):.2f}, max {max(gcs):.2f}\n")
    top_orgs = org_counts.most_common(20)
    summary_lines.append("Top organisms (n):\n")
    for org, n in top_orgs:
        summary_lines.append(f"  - {org}: {n}\n")

    with open(os.path.join(OUT_DIR, 'summary.txt'), 'w') as f:
        f.writelines(summary_lines)

    # Per-sequence table
    write_tsv(os.path.join(OUT_DIR, 'per_sequence_stats.tsv'), per_seq_rows,
              header=['header', 'organism', 'nt_len', 'aa_len', 'gc_percent'])

    # Organism counts
    write_tsv(os.path.join(OUT_DIR, 'organism_counts.tsv'), [(o, n) for o, n in top_orgs], header=['organism', 'count'])

    # Codon usage (aggregate)
    codon_rows = []
    total_codons = sum(agg_codon[c] for c in CODONS)
    for c in CODONS:
        cnt = agg_codon.get(c, 0)
        freq = cnt / total_codons if total_codons else 0.0
        codon_rows.append((c, cnt, f"{freq:.6f}"))
    write_tsv(os.path.join(OUT_DIR, 'codon_usage_aggregate.tsv'), codon_rows, header=['codon', 'count', 'frequency'])

if __name__ == '__main__':
    main()
