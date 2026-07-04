"""
Tier 3 Filter: ESMFold pLDDT
Filters sequences by predicted structural confidence using ESMFold.
Uses original ESMFold from shared environment.
- pLDDT > 70: good structural confidence
- pLDDT > 80: high structural confidence
Runs in hours on ~4K sequences (GPU recommended).
"""
import argparse
import torch
import numpy as np
from pathlib import Path
from Bio import SeqIO
import esm


def get_esmfold_plddt(sequence, model):
    """Run ESMFold on a single sequence and return mean pLDDT."""
    with torch.no_grad():
        # Use the working API from your environment
        output = model.infer_pdb(sequence)
    
    # Extract pLDDT from the PDB string (stored in B-factor column)
    plddt_scores = []
    for line in output.split('\n'):
        if line.startswith('ATOM'):
            # B-factor is in columns 61-66 of PDB format
            try:
                bfactor = float(line[60:66].strip())
                if bfactor > 0:
                    plddt_scores.append(bfactor)
            except (ValueError, IndexError):
                continue
    
    if plddt_scores:
        return float(np.mean(plddt_scores))
    else:
        return 0.0


def filter_by_plddt(input_fasta, output_fasta, min_plddt=70.0, scores_file=None):
    print("Loading ESMFold 150M model (smallest)...")
    # Use the smallest ESMFold model to minimize memory usage
    model = esm.pretrained.esmfold_structure_module_only_150M()
    model = model.eval()
    # Use full precision to avoid LayerNorm errors
    model = model.float()
    print("ESMFold 150M model loaded.")

    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    passed, failed = [], []
    scores = {}

    for i, record in enumerate(sequences):
        seq = str(record.seq)
        if i % 100 == 0:
            print(f"  Processing {i+1}/{len(sequences)}...")
        try:
            plddt = get_esmfold_plddt(seq, model)
            scores[record.id] = plddt
            if plddt >= min_plddt:
                passed.append(record)
            else:
                failed.append(record)
        except Exception as e:
            print(f"  Warning: ESMFold failed for {record.id}: {e}")
            failed.append(record)

    SeqIO.write(passed, output_fasta, "fasta")

    if scores_file:
        with open(scores_file, "w") as f:
            f.write("sequence_id,plddt\n")
            for seq_id, plddt in scores.items():
                f.write(f"{seq_id},{plddt:.2f}\n")
        print(f"  pLDDT scores saved to {scores_file}")

    total = len(sequences)
    print(f"ESMFold pLDDT filter (>{min_plddt}):")
    print(f"  Input:  {total} sequences")
    print(f"  Passed: {len(passed)} ({100*len(passed)/total:.1f}%)")
    print(f"  Failed: {len(failed)} sequences")
    print(f"  Output: {output_fasta}")
    if scores:
        print(f"  Mean pLDDT (passed): {np.mean([scores[r.id] for r in passed if r.id in scores]):.1f}")
    return passed


def main():
    parser = argparse.ArgumentParser(description="Filter sequences by ESMFold pLDDT score.")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    parser.add_argument("--min_plddt", type=float, default=70.0, help="Minimum mean pLDDT (default: 70)")
    parser.add_argument("--scores_file", default=None, help="Optional: save pLDDT scores to CSV")
    args = parser.parse_args()

    filter_by_plddt(args.input, args.output, args.min_plddt, args.scores_file)


if __name__ == "__main__":
    main()
