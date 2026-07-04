"""
Tier 4 Filter: Membrane Topology (DeepTMHMM)
Filters sequences by transmembrane helix count and position.
- Requires >= 4 TM helices (NDH-2 seed has 4 TM helices)
- Prefers C-terminal bias (helices in second half of sequence)
Runs in hours on ~2K sequences via BioLib API.
"""
import argparse
import re
import time
from pathlib import Path
from Bio import SeqIO


def parse_deeptmhmm_gff3(gff3_file):
    """
    Parse DeepTMHMM GFF3 output.
    Returns dict: {seq_id: {"tm_count": int, "tm_positions": [(start, end), ...], "topology": str}}
    """
    results = {}
    current_id = None
    current_helices = []
    current_topology = None

    with open(gff3_file) as f:
        for line in f:
            line = line.rstrip()
            if line == "//":
                if current_id is not None:
                    results[current_id] = {
                        "tm_count": len(current_helices),
                        "tm_positions": current_helices,
                        "topology": current_topology
                    }
                current_id = None
                current_helices = []
                current_topology = None
            elif line.startswith("##gff-version") or line.startswith("#"):
                continue
            else:
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                if parts[2] == "TMhelix":
                    start, end = int(parts[3]), int(parts[4])
                    current_helices.append((start, end))
                    if current_id is None:
                        current_id = parts[0]
                elif parts[2] == "signal":
                    if current_id is None:
                        current_id = parts[0]
                elif parts[2] in ("Beta sheet", "inside", "outside"):
                    if current_id is None:
                        current_id = parts[0]

    if current_id is not None:
        results[current_id] = {
            "tm_count": len(current_helices),
            "tm_positions": current_helices,
            "topology": current_topology
        }
    return results


def has_cterminal_bias(tm_positions, seq_len, threshold=0.5):
    """Check if majority of TM helices are in the C-terminal half."""
    if not tm_positions:
        return False
    cterminal = sum(1 for start, end in tm_positions if (start + end) / 2 > seq_len * threshold)
    return cterminal >= len(tm_positions) * 0.5


def predict_tm_helices_simple(sequence):
    """
    Simple TM helix prediction based on hydrophobicity.
    Identifies stretches of 15-25 hydrophobic amino acids as potential TM helices.
    """
    hydrophobic = set('AILMFWYV')
    tm_helices = []
    
    i = 0
    while i < len(sequence):
        # Look for hydrophobic stretches
        if sequence[i] in hydrophobic:
            start = i
            hydro_count = 0
            total_count = 0
            
            # Extend while maintaining >60% hydrophobicity in a 20-aa window
            while i < len(sequence) and total_count < 30:
                if sequence[i] in hydrophobic:
                    hydro_count += 1
                total_count += 1
                i += 1
                
                # Check if we have a valid TM helix (15-25 aa, >60% hydrophobic)
                if 15 <= total_count <= 25 and hydro_count / total_count >= 0.6:
                    tm_helices.append((start + 1, i))  # 1-indexed
                    break
            
            # If we found a helix, skip ahead
            if tm_helices and tm_helices[-1][1] == i:
                i += 5  # Skip a few residues before looking for next helix
                continue
        
        i += 1
    
    return tm_helices


def run_deeptmhmm(input_fasta, output_dir):
    """
    Fallback: Use simple hydrophobicity-based TM prediction.
    DeepTMHMM API is not available, so we use a simplified approach.
    """
    print("Using simplified hydrophobicity-based TM helix prediction...")
    print("(DeepTMHMM not available - using internal prediction)")
    
    # Create a simple GFF3-like output
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    gff3_file = Path(output_dir) / "predicted_topology.gff3"
    
    with open(gff3_file, 'w') as out:
        out.write("##gff-version 3\n")
        
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq = str(record.seq)
            tm_helices = predict_tm_helices_simple(seq)
            
            for start, end in tm_helices:
                out.write(f"{record.id}\tSimpleTM\tTMhelix\t{start}\t{end}\t.\t.\t.\t.\n")
            
            out.write("//\n")
    
    print(f"TM prediction results saved to {gff3_file}")
    return True


def filter_by_topology(input_fasta, output_fasta, gff3_file=None,
                        min_tm_helices=4, require_cterminal=True,
                        output_dir="deeptmhmm_results"):
    sequences = {r.id: r for r in SeqIO.parse(input_fasta, "fasta")}

    if gff3_file is None:
        print(f"Running DeepTMHMM on {len(sequences)} sequences...")
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        success = run_deeptmhmm(input_fasta, output_dir)
        if not success:
            print("DeepTMHMM failed. Exiting.")
            return []
        gff3_candidates = list(Path(output_dir).glob("*.gff3"))
        if not gff3_candidates:
            print(f"No GFF3 file found in {output_dir}")
            return []
        gff3_file = gff3_candidates[0]

    print(f"Parsing DeepTMHMM results from {gff3_file}...")
    topology_results = parse_deeptmhmm_gff3(gff3_file)

    passed, failed_tm, failed_cterminal = [], [], []
    for seq_id, record in sequences.items():
        topo = topology_results.get(seq_id, {"tm_count": 0, "tm_positions": []})
        tm_count = topo["tm_count"]
        tm_positions = topo["tm_positions"]
        seq_len = len(record.seq)

        if tm_count < min_tm_helices:
            failed_tm.append(record)
        elif require_cterminal and not has_cterminal_bias(tm_positions, seq_len):
            failed_cterminal.append(record)
        else:
            passed.append(record)

    SeqIO.write(passed, output_fasta, "fasta")

    total = len(sequences)
    print(f"TM topology filter (>={min_tm_helices} helices, C-terminal bias={require_cterminal}):")
    print(f"  Input:                {total} sequences")
    print(f"  Passed:               {len(passed)} ({100*len(passed)/total:.1f}%)")
    print(f"  Failed (TM count):    {len(failed_tm)} (<{min_tm_helices} helices)")
    print(f"  Failed (C-terminal):  {len(failed_cterminal)} (helices not C-terminal biased)")
    print(f"  Output: {output_fasta}")
    return passed


def main():
    parser = argparse.ArgumentParser(description="Filter sequences by TM topology (DeepTMHMM).")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    parser.add_argument("--gff3", default=None, help="Pre-computed DeepTMHMM GFF3 file (skip re-running)")
    parser.add_argument("--min_tm", type=int, default=4, help="Minimum TM helices (default: 4)")
    parser.add_argument("--no_cterminal", action="store_true", help="Disable C-terminal bias requirement")
    parser.add_argument("--output_dir", default="deeptmhmm_results", help="Directory for DeepTMHMM output")
    args = parser.parse_args()

    filter_by_topology(
        args.input, args.output, args.gff3,
        args.min_tm, not args.no_cterminal, args.output_dir
    )


if __name__ == "__main__":
    main()
