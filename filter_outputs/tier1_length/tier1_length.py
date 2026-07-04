"""
Tier 1 Filter: Length
Filters sequences by amino acid length (350-800 aa for NDH-2).
Fast: runs in seconds on 10K sequences.
"""
import argparse
from pathlib import Path
from Bio import SeqIO


def filter_by_length(input_fasta, output_fasta, min_len=350, max_len=800):
    passed, failed = [], []
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_len = len(record.seq)
        if min_len <= seq_len <= max_len:
            passed.append(record)
        else:
            failed.append(record)

    SeqIO.write(passed, output_fasta, "fasta")

    print(f"Length filter ({min_len}-{max_len} aa):")
    print(f"  Input:  {len(passed) + len(failed)} sequences")
    print(f"  Passed: {len(passed)} sequences ({100*len(passed)/(len(passed)+len(failed)):.1f}%)")
    print(f"  Failed: {len(failed)} sequences")
    print(f"  Output: {output_fasta}")
    return passed


def main():
    parser = argparse.ArgumentParser(description="Filter sequences by length.")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    parser.add_argument("--min_len", type=int, default=350, help="Minimum length (default: 350)")
    parser.add_argument("--max_len", type=int, default=800, help="Maximum length (default: 800)")
    args = parser.parse_args()

    filter_by_length(args.input, args.output, args.min_len, args.max_len)


if __name__ == "__main__":
    main()
