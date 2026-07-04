#!/usr/bin/env python3
import os
import sys
import re
import math
from collections import Counter, defaultdict
from datetime import datetime

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except Exception:
    matplotlib = None
    plt = None


def ts():
    return datetime.now().isoformat(timespec='seconds')


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


def translate(dna):
    table = {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    }
    aa = []
    for i in range(0, len(dna)-2, 3):
        aa.append(table.get(dna[i:i+3], 'X'))
    return ''.join(aa)


def parse_field(header, key):
    # Extract [key=...] from FASTA header if present
    m = re.search(rf'\[{re.escape(key)}=([^\]]+)\]', header)
    return m.group(1).strip() if m else ''


def gc_pct(seq):
    if not seq:
        return 0.0
    g = seq.count('G') + seq.count('C')
    return 100.0 * g / len(seq)


def kmerset(seq, k=9):
    return {seq[i:i+k] for i in range(0, len(seq)-k+1)} if len(seq) >= k else set()


def jaccard(a, b):
    if not a or not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def main():
    # Args: [input_dir] [parent_cds_fasta] [parent_protein_fasta]
    in_dir = sys.argv[1] if len(sys.argv) > 1 else 'acquire_homologs/outputs_strict_hmm'
    parent_cds_path = sys.argv[2] if len(sys.argv) > 2 else ''
    parent_prot_path = sys.argv[3] if len(sys.argv) > 3 else ''

    cds_fa = os.path.join(in_dir, 'ndh2_cds_expanded.fasta')
    meta_tsv = os.path.join(in_dir, 'ndh2_cds_expanded_metadata.tsv')
    prot_fa = os.path.join(in_dir, 'ndh2_proteins_expanded.fasta')
    plots_dir = os.path.join(in_dir, 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    print(f'[{ts()}] Loading CDS from {cds_fa}')
    entries = list(fasta_stream(cds_fa))
    print(f'[{ts()}] Loaded {len(entries)} sequences')

    # Collect metrics
    nt_lengths = []
    aa_lengths = []
    gcs = []
    organisms = Counter()
    per_seq = []  # (name, nt_len, aa_len, gc, org, jaccard9_dna, jaccard3_protein)

    # Load parent CDS if provided
    parent_seq = ''
    parent_jset = set()
    if parent_cds_path and os.path.exists(parent_cds_path):
        for h, s in fasta_stream(parent_cds_path):
            parent_seq = s
            break
        if parent_seq:
            parent_jset = kmerset(parent_seq, k=9)
            print(f'[{ts()}] Loaded parent CDS: {len(parent_seq)} nt (k-mer set size {len(parent_jset)})')
        else:
            print(f'[{ts()}] WARNING: parent CDS file had no sequence: {parent_cds_path}')
    else:
        if parent_cds_path:
            print(f'[{ts()}] WARNING: parent CDS path not found: {parent_cds_path}')

    # Load parent protein if provided
    parent_prot = ''
    parent_prot_kset = set()
    if parent_prot_path and os.path.exists(parent_prot_path):
        for h, s in fasta_stream(parent_prot_path):
            parent_prot = s
            break
        if parent_prot:
            parent_prot_kset = kmerset(parent_prot, k=3)
            print(f'[{ts()}] Loaded parent protein: {len(parent_prot)} aa (k-mer set size {len(parent_prot_kset)})')
        else:
            print(f'[{ts()}] WARNING: parent protein file had no sequence: {parent_prot_path}')
    else:
        if parent_prot_path:
            print(f'[{ts()}] WARNING: parent protein path not found: {parent_prot_path}')

    for h, s in entries:
        nt = len(s)
        aa = nt // 3
        g = gc_pct(s)
        org = parse_field(h, 'organism') or parse_field(h, 'organism_name')
        nt_lengths.append(nt)
        aa_lengths.append(aa)
        gcs.append(g)
        if org:
            organisms[org] += 1
        jac = jaccard(kmerset(s, 9), parent_jset) if parent_jset else float('nan')
        # protein-level 3-mer Jaccard vs parent protein
        pseq = translate(s)
        pjac = jaccard(kmerset(pseq, 3), parent_prot_kset) if parent_prot_kset else float('nan')
        per_seq.append((h, nt, aa, g, org, jac, pjac))

    # Write per-sequence metrics TSV
    with open(os.path.join(in_dir, 'per_sequence_metrics.tsv'), 'w') as f:
        f.write('header\tnt_len\taa_len\tgc_pct\torganism\tjaccard_k9_to_parent_DNA\tjaccard_k3_to_parent_protein\n')
        for h, nt, aa, g, org, jac, pjac in per_seq:
            f.write(f'{h}\t{nt}\t{aa}\t{g:.3f}\t{org}\t{("%.5f" % jac) if (jac==jac) else ""}\t{("%.5f" % pjac) if (pjac==pjac) else ""}\n')

    if matplotlib is None:
        print(f'[{ts()}] matplotlib not available; wrote metrics TSV, skipping plots')
        return

    # Plots
    # 1) AA length histogram
    plt.figure(figsize=(6,4))
    plt.hist(aa_lengths, bins=40, color='#4C78A8', edgecolor='black', alpha=0.8)
    plt.xlabel('AA length')
    plt.ylabel('Count')
    plt.title('NDH2 AA length distribution')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'aa_length_hist.png'), dpi=200)
    plt.close()

    # 2) GC% histogram
    plt.figure(figsize=(6,4))
    plt.hist(gcs, bins=40, color='#F58518', edgecolor='black', alpha=0.8)
    plt.xlabel('GC%')
    plt.ylabel('Count')
    plt.title('NDH2 GC% distribution')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'gc_pct_hist.png'), dpi=200)
    plt.close()

    # 3) Length vs GC scatter (subsample if huge)
    xs = aa_lengths
    ys = gcs
    if len(xs) > 30000:
        step = max(1, len(xs)//30000)
        xs = xs[::step]
        ys = ys[::step]
    plt.figure(figsize=(6,4))
    plt.scatter(xs, ys, s=4, alpha=0.3, color='#54A24B')
    plt.xlabel('AA length')
    plt.ylabel('GC%')
    plt.title('Length vs GC%')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'length_vs_gc_scatter.png'), dpi=200)
    plt.close()

    # 4) Top organisms bar (top 20)
    if organisms:
        top = organisms.most_common(20)
        labels = [x[0] for x in top]
        counts = [x[1] for x in top]
        plt.figure(figsize=(8,6))
        y_pos = list(range(len(labels)))
        plt.barh(y_pos, counts, color='#E45756', alpha=0.8)
        plt.yticks(y_pos, labels, fontsize=8)
        plt.xlabel('Count')
        plt.title('Top 20 organisms')
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'top_organisms_bar.png'), dpi=200)
        plt.close()

    # 5) Jaccard k-mer DNA similarity histogram (if parent CDS provided)
    if parent_jset:
        vals = [row[5] for row in per_seq if (row[5]==row[5])]
        plt.figure(figsize=(6,4))
        plt.hist(vals, bins=40, color='#72B7B2', edgecolor='black', alpha=0.8)
        plt.xlabel('Jaccard similarity to parent CDS (k=9)')
        plt.ylabel('Count')
        plt.title('DNA similarity to parent (k-mer Jaccard)')
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'dna_similarity_jaccard_k9.png'), dpi=200)
        plt.close()

    # 6) Jaccard k-mer protein similarity histogram (if parent protein provided)
    if parent_prot_kset:
        vals = [row[6] for row in per_seq if (row[6]==row[6])]
        plt.figure(figsize=(6,4))
        plt.hist(vals, bins=40, color='#B279A2', edgecolor='black', alpha=0.8)
        plt.xlabel('Jaccard similarity to parent protein (k=3)')
        plt.ylabel('Count')
        plt.title('Protein similarity to parent (k-mer Jaccard)')
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'protein_similarity_jaccard_k3.png'), dpi=200)
        plt.close()

    print(f'[{ts()}] Plots written to {plots_dir}/')

if __name__ == '__main__':
    main()
