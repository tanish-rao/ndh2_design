#!/usr/bin/env python3
"""
Tier 4 Filter: Simple Membrane Topology Prediction
Filters sequences based on predicted transmembrane helices using hydrophobicity analysis.
This is a fallback when DeepTMHMM is not available.
NDH-2 typically has 4-6 transmembrane helices.
"""

import argparse
import numpy as np
from pathlib import Path
from Bio import SeqIO
import pandas as pd

# Kyte-Doolittle hydropathy index
HYDROPATHY = {
    'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5,
    'M': 1.9, 'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8,
    'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'E': -3.5,
    'Q': -3.5, 'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5
}

def calculate_hydropathy_profile(sequence, window_size=19):
    """Calculate sliding window hydropathy profile."""
    profile = []
    half_window = window_size // 2
    
    for i in range(len(sequence)):
        start = max(0, i - half_window)
        end = min(len(sequence), i + half_window + 1)
        
        window_sum = 0
        for j in range(start, end):
            aa = sequence[j]
            window_sum += HYDROPATHY.get(aa, 0)
        
        profile.append(window_sum / (end - start))
    
    return profile

def predict_tm_helices(sequence, min_length=17, max_length=25, threshold=1.6):
    """Predict transmembrane helices based on hydropathy."""
    profile = calculate_hydropathy_profile(sequence)
    
    tm_helices = []
    in_helix = False
    helix_start = None
    
    for i, score in enumerate(profile):
        if score >= threshold and not in_helix:
            # Start of potential helix
            in_helix = True
            helix_start = i
        elif score < threshold and in_helix:
            # End of helix
            helix_end = i
            helix_length = helix_end - helix_start
            
            if min_length <= helix_length <= max_length:
                tm_helices.append((helix_start, helix_end))
            
            in_helix = False
            helix_start = None
    
    # Check if sequence ends with a helix
    if in_helix and helix_start is not None:
        helix_length = len(sequence) - helix_start
        if min_length <= helix_length <= max_length:
            tm_helices.append((helix_start, len(sequence)))
    
    return tm_helices

def analyze_topology(sequence):
    """Analyze topology of a sequence."""
    tm_helices = predict_tm_helices(sequence)
    
    # Calculate additional metrics
    seq_length = len(sequence)
    tm_coverage = sum(end - start for start, end in tm_helices) / seq_length
    
    # Check for C-terminal bias (more helices in second half)
    first_half = sum(1 for start, end in tm_helices if (start + end) / 2 < seq_length / 2)
    second_half = len(tm_helices) - first_half
    c_terminal_bias = second_half > first_half
    
    return {
        'tm_count': len(tm_helices),
        'tm_positions': tm_helices,
        'tm_coverage': tm_coverage,
        'c_terminal_bias': c_terminal_bias,
        'length': seq_length
    }

def filter_by_topology(input_fasta, output_fasta, min_tm_helices=4, max_tm_helices=8,
                       require_c_bias=False, topology_file=None):
    """Filter sequences based on membrane topology."""
    print(f"Loading sequences from {input_fasta}...")
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    
    print(f"Analyzing topology for {len(sequences)} sequences...")
    
    passed, failed = [], []
    topology_results = []
    
    for i, record in enumerate(sequences):
        seq = str(record.seq)
        if i % 100 == 0:
            print(f"  Processing {i+1}/{len(sequences)}...")
        
        topology = analyze_topology(seq)
        topology['sequence_id'] = record.id
        topology_results.append(topology)
        
        # Apply filters
        tm_count = topology['tm_count']
        passes_tm_filter = min_tm_helices <= tm_count <= max_tm_helices
        passes_bias_filter = not require_c_bias or topology['c_terminal_bias']
        
        if passes_tm_filter and passes_bias_filter:
            passed.append(record)
        else:
            failed.append(record)
    
    # Write results
    SeqIO.write(passed, output_fasta, "fasta")
    
    # Save topology predictions
    if topology_file:
        df = pd.DataFrame(topology_results)
        df.to_csv(topology_file, index=False)
        print(f"  Topology predictions saved to {topology_file}")
    
    # Report results
    total = len(sequences)
    print(f"\nTopology filter ({min_tm_helices}-{max_tm_helices} TM helices")
    if require_c_bias:
        print(f", C-terminal bias required")
    print(f"):")
    print(f"  Input:  {total} sequences")
    print(f"  Passed: {len(passed)} ({100*len(passed)/total:.1f}%)")
    print(f"  Failed: {len(failed)} sequences")
    print(f"  Output: {output_fasta}")
    
    # TM helix distribution
    tm_counts = [r['tm_count'] for r in topology_results]
    print(f"\nTM helix distribution:")
    for tm in sorted(set(tm_counts)):
        count = tm_counts.count(tm)
        passed_count = sum(1 for r in topology_results if r['tm_count'] == tm and r['sequence_id'] in [rec.id for rec in passed])
        print(f"  {tm} TM helices: {count} total ({passed_count} passed)")
    
    # Coverage statistics
    coverages = [r['tm_coverage'] for r in topology_results]
    print(f"\nTM coverage statistics:")
    print(f"  Mean: {np.mean(coverages):.3f}")
    print(f"  Std:  {np.std(coverages):.3f}")
    print(f"  Min:  {np.min(coverages):.3f}")
    print(f"  Max:  {np.max(coverages):.3f}")
    
    # C-terminal bias
    biased_count = sum(1 for r in topology_results if r['c_terminal_bias'])
    print(f"\nC-terminal bias: {biased_count}/{total} ({100*biased_count/total:.1f}%)")
    
    return passed

def main():
    parser = argparse.ArgumentParser(description="Filter sequences by membrane topology")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    parser.add_argument("--min_tm_helices", type=int, default=4,
                       help="Minimum number of TM helices (default: 4)")
    parser.add_argument("--max_tm_helices", type=int, default=8,
                       help="Maximum number of TM helices (default: 8)")
    parser.add_argument("--require_c_bias", action="store_true",
                       help="Require C-terminal bias in TM helices")
    parser.add_argument("--topology_file", default=None,
                       help="Optional: save topology predictions to CSV")
    args = parser.parse_args()

    filter_by_topology(args.input, args.output, args.min_tm_helices,
                      args.max_tm_helices, args.require_c_bias, args.topology_file)

if __name__ == "__main__":
    main()
