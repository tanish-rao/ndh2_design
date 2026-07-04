#!/usr/bin/env python3
"""
Generate structure and run docking for parent enzyme
Compare to Tier 6 candidates
"""
import os
import sys
import torch
from Bio import SeqIO
import esm

def generate_parent_structure():
    """Generate ESMFold structure for parent enzyme"""
    print("=" * 70)
    print("PARENT ENZYME STRUCTURE GENERATION")
    print("=" * 70)
    
    # Load ESMFold model
    print("\nLoading ESMFold model...")
    model = esm.pretrained.esmfold_v1()
    model = model.eval()
    
    # Move to GPU if available
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    model = model.to(device)
    
    # Load parent sequence (resolve relative to this script's directory)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    fasta_file = os.path.join(script_dir, "parent_enzyme.fasta")
    record = list(SeqIO.parse(fasta_file, "fasta"))[0]
    seq_id = record.id
    sequence = str(record.seq)
    
    print(f"\nSequence ID: {seq_id}")
    print(f"Length: {len(sequence)} aa")
    
    # Generate structure
    print("\nGenerating structure with ESMFold...")
    with torch.no_grad():
        output = model.infer_pdb(sequence)
    
    # Save PDB
    output_dir = "pipeline_results/parent_analysis"
    os.makedirs(output_dir, exist_ok=True)
    
    pdb_file = os.path.join(output_dir, f"{seq_id}.pdb")
    with open(pdb_file, "w") as f:
        f.write(output)
    
    print(f"✓ Structure saved to: {pdb_file}")
    
    # Extract pLDDT score
    plddt_scores = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    plddt = float(line[60:66])
                    plddt_scores.append(plddt)
                except (ValueError, IndexError):
                    continue
    
    if plddt_scores:
        mean_plddt = sum(plddt_scores) / len(plddt_scores)
        print(f"✓ Mean pLDDT: {mean_plddt:.2f}")
    
    return pdb_file, mean_plddt

if __name__ == "__main__":
    generate_parent_structure()
    print("\n" + "=" * 70)
    print("STRUCTURE GENERATION COMPLETE")
    print("=" * 70)
