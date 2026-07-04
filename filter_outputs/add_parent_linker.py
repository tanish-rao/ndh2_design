#!/usr/bin/env python3
import argparse
from pathlib import Path
from Bio import SeqIO


def add_linker_to_fasta(input_fasta: Path, parent_fasta: Path, start_aa: int, output_fasta: Path) -> int:
    parent_records = list(SeqIO.parse(str(parent_fasta), 'fasta'))
    if not parent_records:
        raise RuntimeError(f'No records found in parent FASTA: {parent_fasta}')
    parent_seq = str(parent_records[0].seq)
    if start_aa < 1 or start_aa > len(parent_seq):
        raise ValueError(f'start_aa {start_aa} out of range for parent length {len(parent_seq)}')
    tail = parent_seq[start_aa-1:]

    input_records = list(SeqIO.parse(str(input_fasta), 'fasta'))
    out_records = []
    for rec in input_records:
        new_seq = str(rec.seq) + tail
        rec_out = rec[:]
        rec_out.seq = rec.seq.__class__(new_seq)
        out_records.append(rec_out)

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(out_records, str(output_fasta), 'fasta')
    return len(out_records)


def main():
    ap = argparse.ArgumentParser(description='Append parent linker domain to sequences')
    ap.add_argument('--input_fasta', required=True, help='Input protein FASTA (e.g., Tier 2 output)')
    ap.add_argument('--parent_fasta', required=True, help='Parent enzyme FASTA')
    ap.add_argument('--start_aa', type=int, default=423, help='1-indexed AA position to take tail from parent (default 423)')
    ap.add_argument('--output_fasta', required=True, help='Output protein FASTA with linker appended')
    args = ap.parse_args()

    n = add_linker_to_fasta(Path(args.input_fasta), Path(args.parent_fasta), args.start_aa, Path(args.output_fasta))
    print(f'Wrote {n} sequences with appended linker to {args.output_fasta}')


if __name__ == '__main__':
    main()
