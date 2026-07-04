#!/usr/bin/env python3
import argparse
import os
import json
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def load_scores(csv_path: Path) -> Dict[str, float]:
    import csv
    scores: Dict[str, float] = {}
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        # Expect columns: sequence_id,plddt
        # Fallbacks: id, score
        id_key = 'sequence_id'
        score_key = 'plddt'
        header = reader.fieldnames or []
        if 'sequence_id' not in header and 'id' in header:
            id_key = 'id'
        if 'plddt' not in header and 'score' in header:
            score_key = 'score'
        for row in reader:
            sid = row.get(id_key)
            val = row.get(score_key)
            if sid is None or val is None:
                continue
            try:
                scores[sid] = float(val)
            except ValueError:
                continue
    return scores


def top_n(scores: Dict[str, float], n: int) -> List[Tuple[str, float]]:
    return sorted(scores.items(), key=lambda kv: kv[1], reverse=True)[:n]


def main():
    ap = argparse.ArgumentParser(description='Compare pre-linker vs linker-appended Tier 3 pLDDT')
    ap.add_argument('--pre_csv', required=True, help='Pre-linker pLDDT CSV (tier3_plddt_scores.csv)')
    ap.add_argument('--post_csv', required=True, help='Post-linker pLDDT CSV (tier3_plddt_scores.csv)')
    ap.add_argument('--out_dir', required=True, help='Directory to save comparison outputs')
    ap.add_argument('--top_n', type=int, default=100, help='Top-N to compare (default 100)')
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pre = load_scores(Path(args.pre_csv))
    post = load_scores(Path(args.post_csv))

    # Overlay histograms
    pre_vals = list(pre.values())
    post_vals = list(post.values())
    if pre_vals and post_vals:
        plt.figure(figsize=(7,5))
        plt.hist(pre_vals, bins=40, alpha=0.5, density=True, label='Pre-linker', color='#4e79a7')
        plt.hist(post_vals, bins=40, alpha=0.5, density=True, label='Linker-appended', color='#e15759')
        plt.xlabel('pLDDT')
        plt.ylabel('Density')
        plt.title('Tier 3 pLDDT distribution: pre-linker vs linker-appended')
        plt.legend()
        plt.tight_layout()
        plt.savefig(out_dir / 'plddt_hist_overlay.png', dpi=150)
        plt.close()

    # Scatter for intersecting IDs
    common = sorted(set(pre.keys()) & set(post.keys()))
    x = [pre[sid] for sid in common]
    y = [post[sid] for sid in common]
    if common:
        import numpy as np
        from scipy.stats import pearsonr, spearmanr
        plt.figure(figsize=(6,6))
        plt.scatter(x, y, s=8, alpha=0.4, color='#59a14f')
        lims = [min(x+y)-1, max(x+y)+1]
        plt.plot(lims, lims, 'k--', linewidth=1)
        plt.xlabel('Pre-linker pLDDT')
        plt.ylabel('Linker-appended pLDDT')
        plt.title('Tier 3 pLDDT: per-sequence comparison')
        plt.tight_layout()
        plt.savefig(out_dir / 'plddt_scatter_pre_vs_post.png', dpi=150)
        plt.close()
        try:
            pcc = pearsonr(x, y)[0]
            scc = spearmanr(x, y)[0]
        except Exception:
            pcc = None
            scc = None
    else:
        pcc = None
        scc = None

    # Compare top-N overlap
    tn = int(args.top_n)
    pre_top = top_n(pre, tn)
    post_top = top_n(post, tn)
    pre_ids = [sid for sid, _ in pre_top]
    post_ids = [sid for sid, _ in post_top]
    overlap = sorted(set(pre_ids) & set(post_ids))

    # Also check top-1 IDs
    top1_pre = pre_top[0][0] if pre_top else None
    top1_post = post_top[0][0] if post_top else None

    # Save summary JSON
    summary = {
        'pre_count': len(pre),
        'post_count': len(post),
        'common_ids': len(common),
        'pearson_r': pcc,
        'spearman_r': scc,
        'top_n': tn,
        'top1_pre': top1_pre,
        'top1_post': top1_post,
        'top_overlap_count': len(overlap),
        'top_overlap_jaccard': (len(overlap) / (2*tn - len(overlap))) if tn > 0 else None,
    }
    with open(out_dir / 'comparison_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)

    # Save overlap lists
    with open(out_dir / 'top_pre_ids.txt', 'w') as f:
        for sid in pre_ids:
            f.write(sid + '\n')
    with open(out_dir / 'top_post_ids.txt', 'w') as f:
        for sid in post_ids:
            f.write(sid + '\n')
    with open(out_dir / 'top_overlap_ids.txt', 'w') as f:
        for sid in overlap:
            f.write(sid + '\n')

    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
