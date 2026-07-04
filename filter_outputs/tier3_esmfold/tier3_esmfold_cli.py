#!/usr/bin/env python3
"""
Tier 3 Filter: ESMFold pLDDT using CLI
Filters sequences by predicted structural confidence using ESMFold CLI.
- pLDDT > 70: good structural confidence
- pLDDT > 80: high structural confidence
Processes sequences in batches to avoid memory issues.
"""
import argparse
import os
import subprocess
import tempfile
import glob
from pathlib import Path
from Bio import SeqIO
import pandas as pd


def extract_plddt_from_pdb(pdb_file):
    """Extract mean pLDDT from a PDB file."""
    plddt_scores = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                # B-factor column contains pLDDT score
                try:
                    bfactor = float(line[60:66].strip())
                    if bfactor > 0:
                        plddt_scores.append(bfactor)
                except:
                    continue
    
    if plddt_scores:
        return sum(plddt_scores) / len(plddt_scores)
    else:
        return 0.0


def process_batch(sequences, batch_num, output_dir, max_tokens_per_batch=512):
    """Process a batch of sequences using ESMFold CLI."""
    # Create temporary FASTA for this batch
    batch_file = f"batch_{batch_num}.fasta"
    with open(batch_file, 'w') as f:
        for seq_id, seq in sequences:
            f.write(f">{seq_id}\n{seq}\n")
    
    # Run ESMFold CLI
    cmd = [
        "python3", "-m", "esm.scripts.fold",
        "-i", batch_file,
        "-o", output_dir,
        "--max-tokens-per-batch", str(max_tokens_per_batch)
    ]
    
    print(f"  Processing batch {batch_num} ({len(sequences)} sequences)...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"  Error in batch {batch_num}: {result.stderr}")
        return {}
    
    # Clean up batch file
    os.remove(batch_file)
    
    # Extract pLDDT scores from generated PDB files
    scores = {}
    for seq_id, _ in sequences:
        pdb_file = os.path.join(output_dir, f"{seq_id}.pdb")
        if os.path.exists(pdb_file):
            plddt = extract_plddt_from_pdb(pdb_file)
            scores[seq_id] = plddt
    
    return scores


def filter_by_plddt(input_fasta, output_fasta, min_plddt=70.0, scores_file=None, batch_size=10):
    """Filter sequences by ESMFold pLDDT score using CLI."""
    print(f"Loading sequences from {input_fasta}...")
    sequences = [(record.id, str(record.seq)) for record in SeqIO.parse(input_fasta, "fasta")]
    
    print(f"Processing {len(sequences)} sequences in batches of {batch_size}...")
    
    # Create output directory
    output_dir = "pipeline_results/esmfold_temp"
    os.makedirs(output_dir, exist_ok=True)
    
    all_scores = {}
    passed, failed = [], []
    
    # Process in batches
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i+batch_size]
        batch_num = i // batch_size + 1
        
        scores = process_batch(batch, batch_num, output_dir)
        all_scores.update(scores)
        
        # Filter based on pLDDT
        for seq_id, seq in batch:
            if seq_id in scores:
                if scores[seq_id] >= min_plddt:
                    # Find the original record
                    for record in SeqIO.parse(input_fasta, "fasta"):
                        if record.id == seq_id:
                            passed.append(record)
                            break
                else:
                    failed.append(seq_id)
            else:
                failed.append(seq_id)
                print(f"  Warning: No score for {seq_id}")
    
    # Write passed sequences
    SeqIO.write(passed, output_fasta, "fasta")
    
    # Save scores
    if scores_file:
        df = pd.DataFrame(list(all_scores.items()), columns=['sequence_id', 'plddt'])
        df.to_csv(scores_file, index=False)
        print(f"  pLDDT scores saved to {scores_file}")
    
    # Clean up temporary PDB files
    pdb_files = glob.glob(os.path.join(output_dir, "*.pdb"))
    for pdb_file in pdb_files:
        os.remove(pdb_file)
    os.rmdir(output_dir)
    
    # Report results
    total = len(sequences)
    print(f"\nESMFold pLDDT filter (>{min_plddt}):")
    print(f"  Input:  {total} sequences")
    print(f"  Passed: {len(passed)} ({100*len(passed)/total:.1f}%)")
    print(f"  Failed: {len(failed)} sequences")
    print(f"  Output: {output_fasta}")
    
    if all_scores:
        passed_scores = [all_scores[r.id] for r in passed if r.id in all_scores]
        if passed_scores:
            print(f"  Mean pLDDT (passed): {sum(passed_scores)/len(passed_scores):.1f}")
    
    return passed


def main():
    parser = argparse.ArgumentParser(description="Filter sequences by ESMFold pLDDT score using CLI.")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    parser.add_argument("--min_plddt", type=float, default=70.0, help="Minimum mean pLDDT (default: 70)")
    parser.add_argument("--scores_file", default=None, help="Optional: save pLDDT scores to CSV")
    parser.add_argument("--batch_size", type=int, default=10, help="Batch size for processing (default: 10)")
    args = parser.parse_args()

    filter_by_plddt(args.input, args.output, args.min_plddt, args.scores_file, args.batch_size)


if __name__ == "__main__":
    main()
