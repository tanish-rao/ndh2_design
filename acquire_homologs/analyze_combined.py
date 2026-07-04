#!/usr/bin/env python3
import os
import sys
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_meta(path):
    rows=[]
    with open(path,'r') as f:
        header=f.readline().rstrip('\n').split('\t')
        for line in f:
            parts=line.rstrip('\n').split('\t')
            if len(parts)<6:
                continue
            rid=parts[0]
            src=parts[1]
            try: nt=int(parts[2])
            except: nt=math.nan
            try: aa=int(parts[3])
            except: aa=math.nan
            try: gc=float(parts[4]) if parts[4] else math.nan
            except: gc=math.nan
            try: pjac=float(parts[5]) if parts[5] else math.nan
            except: pjac=math.nan
            rows.append((rid,src,nt,aa,gc,pjac))
    return rows


def main():
    # Args: [combined_dir]
    cdir = sys.argv[1] if len(sys.argv)>1 else 'acquire_homologs/combined'
    meta_path=os.path.join(cdir,'combined_metadata.tsv')
    plots=os.path.join(cdir,'plots')
    os.makedirs(plots, exist_ok=True)
    rows=read_meta(meta_path)

    aa=[r[3] for r in rows if not math.isnan(r[3])]
    gc=[r[4] for r in rows if not math.isnan(r[4])]
    pjac=[r[5] for r in rows if not math.isnan(r[5])]

    plt.figure(figsize=(6,4))
    plt.hist(aa, bins=40, color='#4C78A8', edgecolor='black', alpha=0.8)
    plt.xlabel('AA length')
    plt.ylabel('Count')
    plt.title('Combined set: AA length')
    plt.tight_layout()
    plt.savefig(os.path.join(plots,'combined_aa_length_hist.png'), dpi=200)
    plt.close()

    plt.figure(figsize=(6,4))
    plt.hist(gc, bins=40, color='#F58518', edgecolor='black', alpha=0.8)
    plt.xlabel('GC%')
    plt.ylabel('Count')
    plt.title('Combined set: GC%')
    plt.tight_layout()
    plt.savefig(os.path.join(plots,'combined_gc_pct_hist.png'), dpi=200)
    plt.close()

    if pjac:
        plt.figure(figsize=(6,4))
        plt.hist(pjac, bins=40, color='#B279A2', edgecolor='black', alpha=0.8)
        plt.xlabel('Jaccard similarity to parent protein (k=3)')
        plt.ylabel('Count')
        plt.title('Combined set: Protein similarity')
        plt.tight_layout()
        plt.savefig(os.path.join(plots,'combined_protein_similarity_jaccard_k3.png'), dpi=200)
        plt.close()

    # Length vs GC scatter
    if aa and gc:
        plt.figure(figsize=(5.5,4.5))
        plt.scatter(aa, gc, s=6, alpha=0.5, color='#72B7B2', edgecolors='none')
        plt.xlabel('AA length')
        plt.ylabel('GC%')
        plt.title('Combined set: Length vs GC%')
        plt.tight_layout()
        plt.savefig(os.path.join(plots,'combined_length_vs_gc_scatter.png'), dpi=200)
        plt.close()

    print(f"[combined-analysis] N={len(rows)}; Plots -> {plots}")

if __name__=='__main__':
    main()
