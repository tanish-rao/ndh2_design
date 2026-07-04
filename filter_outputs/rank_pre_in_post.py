#!/usr/bin/env python3
import argparse
from pathlib import Path
from typing import Dict, List
from Bio import SeqIO


def load_ids_from_fasta(path: Path) -> List[str]:
    return [rec.id for rec in SeqIO.parse(str(path), 'fasta')]


def load_post_scores(path: Path) -> Dict[str, float]:
    import csv
    scores: Dict[str, float] = {}
    with open(path) as f:
        r = csv.DictReader(f)
        id_key = 'sequence_id' if 'sequence_id' in (r.fieldnames or []) else 'id'
        score_key = 'plddt' if 'plddt' in (r.fieldnames or []) else 'score'
        for row in r:
            sid = row.get(id_key)
            val = row.get(score_key)
            if sid is None or val is None:
                continue
            try:
                scores[sid] = float(val)
            except ValueError:
                continue
    return scores


def main():
    ap = argparse.ArgumentParser(description='Rank pre-linker pass IDs in post-linker ranking')
    ap.add_argument('--pre_pass_fasta', required=True, help='Pre-linker Tier3 passed FASTA (379 seqs)')
    ap.add_argument('--post_scores_csv', required=True, help='Post-linker pLDDT CSV for all seqs')
    args = ap.parse_args()

    pre_ids = load_ids_from_fasta(Path(args.pre_pass_fasta))
    post_scores = load_post_scores(Path(args.post_scores_csv))

    # Build post ranking: sort desc by pLDDT
    ranked = sorted(post_scores.items(), key=lambda kv: kv[1], reverse=True)
    rank_index: Dict[str, int] = {sid: i+1 for i, (sid, _) in enumerate(ranked)}  # 1-based ranks

    # Compute ranks for pre_ids
    ranks = [(sid, rank_index.get(sid, None), post_scores.get(sid, None)) for sid in pre_ids]
    missing = [sid for sid, rk, _ in ranks if rk is None]
    present = [(sid, rk, sc) for sid, rk, sc in ranks if rk is not None]

    if present:
        worst = max(present, key=lambda x: x[1])
        best = min(present, key=lambda x: x[1])
        print({
            'n_pre_pass_ids': len(pre_ids),
            'n_found_in_post': len(present),
            'n_missing_in_post': len(missing),
            'worst_id': worst[0],
            'worst_rank_in_post': worst[1],
            'worst_post_plddt': worst[2],
            'best_id': best[0],
            'best_rank_in_post': best[1],
            'best_post_plddt': best[2],
        })
    else:
        print({'n_pre_pass_ids': len(pre_ids), 'n_found_in_post': 0, 'n_missing_in_post': len(missing)})


if __name__ == '__main__':
    main()
