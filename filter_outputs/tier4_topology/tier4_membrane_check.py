#!/usr/bin/env python3
"""
Tier 4 Filter: Check for membrane-associated features in NDH-2 sequences
L. plantarum NDH-2 should have C-terminal transmembrane helices
"""
import argparse
from pathlib import Path
from Bio import SeqIO
import pandas as pd

def check_hydrophobic_stretches(sequence, min_length=15, threshold=0.6):
    """Check for hydrophobic stretches typical of TM helices."""
    hydrophobic = set('AILMFWYV')
    
    stretches = []
    current_stretch = 0
    
    for aa in sequence:
        if aa in hydrophobic:
            current_stretch += 1
        else:
            if current_stretch >= min_length:
                stretches.append(current_stretch)
            current_stretch = 0
    
    # Check last stretch
    if current_stretch >= min_length:
        stretches.append(current_stretch)
    
    return stretches

def has_cterminal_hydrophobic_region(sequence, last_n_aa=100, min_length=20):
    """Check if the C-terminus has a hydrophobic region."""
    if len(sequence) < last_n_aa:
        return False, 0
    
    c_terminal = sequence[-last_n_aa:]
    stretches = check_hydrophobic_stretches(c_terminal, min_length=min_length)
    
    return len(stretches) > 0, max(stretches) if stretches else 0

def analyze_membrane_features(sequence):
    """Analyze sequence for membrane protein features."""
    length = len(sequence)
    
    # Check overall hydrophobic content
    hydrophobic = set('AILMFWYV')
    hydrophobic_content = sum(1 for aa in sequence if aa in hydrophobic) / length
    
    # Check for long hydrophobic stretches
    all_stretches = check_hydrophobic_stretches(sequence)
    max_stretch = max(all_stretches) if all_stretches else 0
    total_hydrophobic_in_stretches = sum(all_stretches)
    
    # Check C-terminal region specifically
    has_cterm, cterm_length = has_cterminal_hydrophobic_region(sequence)
    
    # Calculate what fraction of sequence is in long hydrophobic stretches
    stretch_fraction = total_hydrophobic_in_stretches / length if length > 0 else 0
    
    return {
        'length': length,
        'hydrophobic_content': hydrophobic_content,
        'max_hydrophobic_stretch': max_stretch,
        'total_hydrophobic_stretches': len(all_stretches),
        'has_cterminal_region': has_cterm,
        'cterminal_max_stretch': cterm_length,
        'stretch_fraction': stretch_fraction
    }

def filter_for_membrane_associated(input_fasta, output_fasta, 
                                   min_hydrophobic_content=0.35,
                                   min_max_stretch=15,
                                   require_cterminal=True,
                                   min_cterminal_stretch=18,
                                   analysis_file=None):
    """Filter for sequences that could be membrane-associated NDH-2."""
    
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    print(f"Analyzing {len(sequences)} sequences for membrane features...")
    
    passed, failed = [], []
    analyses = []
    
    for i, record in enumerate(sequences):
        if i % 100 == 0:
            print(f"  Processing {i+1}/{len(sequences)}...")
        
        seq = str(record.seq)
        features = analyze_membrane_features(seq)
        features['sequence_id'] = record.id
        analyses.append(features)
        
        # Apply filters
        passes_hydrophobic = features['hydrophobic_content'] >= min_hydrophobic_content
        passes_stretch = features['max_hydrophobic_stretch'] >= min_max_stretch
        passes_cterminal = not require_cterminal or features['has_cterminal_region']
        
        if passes_hydrophobic and passes_stretch and passes_cterminal:
            passed.append(record)
        else:
            failed.append(record)
    
    # Write results
    SeqIO.write(passed, output_fasta, "fasta")
    
    # Save analysis
    if analysis_file:
        df = pd.DataFrame(analyses)
        df.to_csv(analysis_file, index=False)
        print(f"  Analysis saved to {analysis_file}")
    
    # Report results
    total = len(sequences)
    print(f"\nMembrane-associated filter results:")
    print(f"  Input: {total} sequences")
    print(f"  Passed: {len(passed)} ({100*len(passed)/total:.1f}%)")
    print(f"  Failed: {len(failed)} sequences")
    print(f"  Output: {output_fasta}")
    
    # Statistics
    df = pd.DataFrame(analyses)
    print(f"\nFeature statistics (all sequences):")
    print(f"  Mean hydrophobic content: {df['hydrophobic_content'].mean():.3f}")
    print(f"  Mean max hydrophobic stretch: {df['max_hydrophobic_stretch'].mean():.1f}")
    print(f"  Sequences with C-terminal region: {df['has_cterminal_region'].sum()} ({100*df['has_cterminal_region'].sum()/total:.1f}%)")
    
    if len(passed) > 0:
        passed_ids = {r.id for r in passed}
        passed_df = df[df['sequence_id'].isin(passed_ids)]
        print(f"\nFeature statistics (passed sequences):")
        print(f"  Mean hydrophobic content: {passed_df['hydrophobic_content'].mean():.3f}")
        print(f"  Mean max hydrophobic stretch: {passed_df['max_hydrophobic_stretch'].mean():.1f}")
        print(f"  Sequences with C-terminal region: {passed_df['has_cterminal_region'].sum()} ({100*passed_df['has_cterminal_region'].sum()/len(passed):.1f}%)")
    
    return passed

def main():
    parser = argparse.ArgumentParser(description="Filter for membrane-associated NDH-2 sequences")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    parser.add_argument("--min_hydrophobic_content", type=float, default=0.35,
                       help="Minimum hydrophobic content (default: 0.35)")
    parser.add_argument("--min_max_stretch", type=int, default=15,
                       help="Minimum longest hydrophobic stretch (default: 15)")
    parser.add_argument("--require_cterminal", action="store_true", default=True,
                       help="Require C-terminal hydrophobic region")
    parser.add_argument("--min_cterminal_stretch", type=int, default=18,
                       help="Minimum C-terminal hydrophobic stretch (default: 18)")
    parser.add_argument("--analysis_file", default=None,
                       help="Save analysis to CSV file")
    args = parser.parse_args()

    filter_for_membrane_associated(args.input, args.output,
                                   args.min_hydrophobic_content,
                                   args.min_max_stretch,
                                   args.require_cterminal,
                                   args.min_cterminal_stretch,
                                   args.analysis_file)

if __name__ == "__main__":
    main()
