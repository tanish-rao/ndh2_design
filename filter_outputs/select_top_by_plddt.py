#!/usr/bin/env python3
import argparse
import json
import re
from pathlib import Path
from typing import Dict, List, Tuple
from Bio import SeqIO

IDX_RE = re.compile(r".*?(\d+)$")

def parse_index(seq_id: str) -> int:
    m = IDX_RE.match(seq_id.replace('>', '').strip())
    if not m:
        raise ValueError(f"Cannot parse numeric index from id: {seq_id}")
    return int(m.group(1))


def load_scores(csv_path: Path) -> List[Tuple[str, float]]:
    import csv
    rows: List[Tuple[str, float]] = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        # Expected columns: sequence_id, plddt
        # Fallbacks if headers differ
        id_key = 'sequence_id' if 'sequence_id' in (reader.fieldnames or []) else 'id'
        score_key = 'plddt' if 'plddt' in (reader.fieldnames or []) else 'score'
        for row in reader:
            sid = row.get(id_key)
            sval = row.get(score_key)
            if sid is None or sval is None:
                continue
            try:
                rows.append((sid, float(sval)))
            except ValueError:
                continue
    # Sort descending by score
    rows.sort(key=lambda kv: kv[1], reverse=True)
    return rows


def select_top(scores: List[Tuple[str, float]], n: int) -> List[Tuple[str, float]]:
    return scores[:n]


def write_subset_proteins(all_fasta: Path, keep_ids: List[str], out_fasta: Path) -> int:
    keep_set = set(keep_ids)
    recs = [r for r in SeqIO.parse(str(all_fasta), 'fasta') if r.id in keep_set]
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(recs, str(out_fasta), 'fasta')
    return len(recs)


def write_subset_codons(source_codons: Path, keep_ids: List[str], out_codons: Path) -> int:
    # Map by numeric index suffix in ID (e.g., gen_872 -> 872)
    idxs = []
    for sid in keep_ids:
        try:
            idxs.append(parse_index(sid))
        except Exception:
            continue
    idx_set = set(idxs)
    lines = [ln.strip() for ln in open(source_codons) if ln.strip()]
    with open(out_codons, 'w') as f:
        kept = 0
        for i, line in enumerate(lines):
            if i in idx_set:
                f.write(line + '\n')
                kept += 1
    return kept


def main():
    ap = argparse.ArgumentParser(description='Select top-N sequences by pLDDT and write proteins/codons subsets')
    ap.add_argument('--scores_csv', required=True, help='CSV with columns including sequence_id and pLDDT')
    ap.add_argument('--proteins_fasta', required=True, help='Protein FASTA to subset (e.g., linker-appended FASTA)')
    ap.add_argument('--source_codons', required=True, help='Original cleaned codon tokens (one seq per line)')
    ap.add_argument('--out_dir', required=True, help='Output directory for subset')
    ap.add_argument('--top_n', type=int, default=500, help='How many top sequences to select (default 500)')
    args = ap.parse_args()

    scores = load_scores(Path(args.scores_csv))
    top = select_top(scores, int(args.top_n))
    top_ids = [sid for sid, _ in top]
    top_scores = [sc for _, sc in top]

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Write proteins subset
    n_prot = write_subset_proteins(Path(args.proteins_fasta), top_ids, out_dir / 'proteins.fasta')
    # Write codons subset
    n_cod = write_subset_codons(Path(args.source_codons), top_ids, out_dir / 'codons.txt')

    summary = {
        'top_n': len(top_ids),
        'min_plddt': min(top_scores) if top_scores else None,
        'max_plddt': max(top_scores) if top_scores else None,
        'proteins_written': n_prot,
        'codons_written': n_cod,
    }
    with open(out_dir / 'summary.json', 'w') as f:
        json.dump(summary, f, indent=2)

    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
