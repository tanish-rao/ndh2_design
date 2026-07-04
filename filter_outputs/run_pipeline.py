"""
NDH-2 Sequence Filtering Pipeline
Tiered filtering pipeline for GenSLM-generated NDH-2 sequences.

Tier 1 (seconds):  Length filter (350-800 aa)
Tier 2 (minutes):  BLAST similarity (30-90% identity to reference)
Tier 3 (hours):    ESMFold pLDDT (>70 structural confidence)
Tier 4 (hours):    TM topology (>=4 helices, C-terminal bias)
Tier 5 (minutes):  Clustering (70% identity, diversity selection)
Tier 6 (days):     Molecular docking (NADH + DHNA binding)
Tier 7 (weeks):    DFT (HOMO-LUMO gap, binding site electronics)

Usage:
  python run_pipeline.py --input generated_ndh2_10k_cpu.fasta --output_dir pipeline_results/
  python run_pipeline.py --input generated_ndh2_10k_cpu.fasta --output_dir pipeline_results/ --tiers 1,2,3
"""
import argparse
import os
import sys
import json
from pathlib import Path
from datetime import datetime

# Add this directory and tier subdirectories to path
base_dir = os.path.dirname(__file__)
sys.path.insert(0, base_dir)
sys.path.insert(0, os.path.join(base_dir, 'tier1_length'))
sys.path.insert(0, os.path.join(base_dir, 'tier2_blast'))
sys.path.insert(0, os.path.join(base_dir, 'tier3_esmfold'))
sys.path.insert(0, os.path.join(base_dir, 'tier4_topology'))
sys.path.insert(0, os.path.join(base_dir, 'tier5_clustering'))
sys.path.insert(0, os.path.join(base_dir, 'tier6_docking'))
sys.path.insert(0, os.path.join(base_dir, 'tier7_dft'))

from tier1_length import filter_by_length
from tier2_blast import filter_by_blast
from tier3_esmfold import filter_by_plddt
from tier4_topology import filter_by_topology
from tier5_clustering import filter_by_clustering
from tier6_docking import filter_by_docking
from tier7_dft import filter_by_dft


def count_seqs(fasta_file):
    """Count sequences in a FASTA file."""
    if not Path(fasta_file).exists():
        return 0
    count = 0
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def log_tier(log, tier, input_file, output_file, start_time, end_time):
    """Log tier results."""
    duration = (end_time - start_time).total_seconds()
    entry = {
        "tier": tier,
        "input_count": count_seqs(input_file),
        "output_count": count_seqs(output_file),
        "duration_seconds": duration,
        "input_file": str(input_file),
        "output_file": str(output_file),
        "timestamp": end_time.isoformat()
    }
    log.append(entry)
    retention = 100 * entry["output_count"] / entry["input_count"] if entry["input_count"] > 0 else 0
    print(f"\n{'='*60}")
    print(f"Tier {tier} complete: {entry['input_count']} → {entry['output_count']} sequences ({retention:.1f}% retained)")
    print(f"Duration: {duration:.1f}s")
    print(f"{'='*60}\n")
    return entry


def main():
    parser = argparse.ArgumentParser(description="NDH-2 tiered filtering pipeline.")
    parser.add_argument("--input", required=True, help="Input FASTA file (generated sequences)")
    parser.add_argument("--output_dir", required=True, help="Output directory for pipeline results")
    parser.add_argument("--tiers", default="1,2,3,4,5,6,7",
                        help="Comma-separated list of tiers to run (default: 1,2,3,4,5,6,7)")

    # Tier 1: Length
    parser.add_argument("--min_len", type=int, default=350, help="Min sequence length (default: 350)")
    parser.add_argument("--max_len", type=int, default=800, help="Max sequence length (default: 800)")

    # Tier 2: BLAST
    parser.add_argument("--ref_fasta", default=None,
                        help="Reference FASTA for BLAST (default: close Firmicutes set)")
    parser.add_argument("--min_identity", type=float, default=30.0, help="Min BLAST identity %% (default: 30)")
    parser.add_argument("--max_identity", type=float, default=90.0, help="Max BLAST identity %% (default: 90)")
    parser.add_argument("--blast_threads", type=int, default=8, help="BLAST threads (default: 8)")

    # Tier 3: ESMFold
    parser.add_argument("--min_plddt", type=float, default=70.0, help="Min ESMFold pLDDT (default: 70)")

    # Tier 4: Topology
    parser.add_argument("--min_tm", type=int, default=4, help="Min TM helices (default: 4)")
    parser.add_argument("--gff3", default=None, help="Pre-computed DeepTMHMM GFF3 file")
    parser.add_argument("--no_cterminal", action="store_true", help="Disable C-terminal bias requirement")

    # Tier 5: Clustering
    parser.add_argument("--cluster_identity", type=float, default=0.7,
                        help="Clustering identity threshold (default: 0.7)")

    # Tier 6: Docking
    parser.add_argument("--structures_dir", default=None, help="Directory with ESMFold PDB structures")
    parser.add_argument("--min_nadh_energy", type=float, default=-7.0, help="Min NADH binding energy kcal/mol")
    parser.add_argument("--min_dhna_energy", type=float, default=-7.0, help="Min DHNA binding energy kcal/mol")

    # Tier 7: DFT
    parser.add_argument("--docking_dir", default=None, help="Directory with docking pose PDB files")
    parser.add_argument("--min_gap", type=float, default=2.0, help="Min HOMO-LUMO gap in eV (default: 2.0)")

    args = parser.parse_args()

    # Parse tiers to run
    tiers_to_run = set(int(t.strip()) for t in args.tiers.split(","))

    # Setup output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Default reference FASTA
    ref_fasta = args.ref_fasta or str(
        Path(__file__).parent.parent /
        "ndh2_ref/model/ref_close_firmicutes_len350-800_dedup.fasta"
    )

    log = []
    current_input = args.input

    print(f"\n{'='*60}")
    print(f"NDH-2 Filtering Pipeline")
    print(f"Input: {args.input} ({count_seqs(args.input)} sequences)")
    print(f"Tiers: {sorted(tiers_to_run)}")
    print(f"Output dir: {output_dir}")
    print(f"{'='*60}\n")

    # -------------------------------------------------------------------------
    # Tier 1: Length
    # -------------------------------------------------------------------------
    if 1 in tiers_to_run:
        print(f"[Tier 1] Length filter ({args.min_len}-{args.max_len} aa)...")
        t1_output = str(output_dir / "tier1_length.fasta")
        start = datetime.now()
        filter_by_length(current_input, t1_output, args.min_len, args.max_len)
        log_tier(log, 1, current_input, t1_output, start, datetime.now())
        current_input = t1_output

    # -------------------------------------------------------------------------
    # Tier 2: BLAST
    # -------------------------------------------------------------------------
    if 2 in tiers_to_run:
        print(f"[Tier 2] BLAST similarity filter ({args.min_identity}-{args.max_identity}% identity)...")
        t2_output = str(output_dir / "tier2_blast.fasta")
        start = datetime.now()
        filter_by_blast(current_input, t2_output, ref_fasta,
                        args.min_identity, args.max_identity, args.blast_threads)
        log_tier(log, 2, current_input, t2_output, start, datetime.now())
        current_input = t2_output

    # -------------------------------------------------------------------------
    # Tier 3: ESMFold pLDDT
    # -------------------------------------------------------------------------
    if 3 in tiers_to_run:
        print(f"[Tier 3] ESMFold pLDDT filter (>{args.min_plddt})...")
        t3_output = str(output_dir / "tier3_esmfold.fasta")
        t3_scores = str(output_dir / "tier3_plddt_scores.csv")
        start = datetime.now()
        filter_by_plddt(current_input, t3_output, args.min_plddt, t3_scores)
        log_tier(log, 3, current_input, t3_output, start, datetime.now())
        current_input = t3_output

    # -------------------------------------------------------------------------
    # Tier 4: TM Topology
    # -------------------------------------------------------------------------
    if 4 in tiers_to_run:
        print(f"[Tier 4] TM topology filter (>={args.min_tm} helices)...")
        t4_output = str(output_dir / "tier4_topology.fasta")
        t4_deeptmhmm_dir = str(output_dir / "tier4_deeptmhmm")
        start = datetime.now()
        filter_by_topology(current_input, t4_output, args.gff3,
                           args.min_tm, not args.no_cterminal, t4_deeptmhmm_dir)
        log_tier(log, 4, current_input, t4_output, start, datetime.now())
        current_input = t4_output

    # -------------------------------------------------------------------------
    # Tier 5: Clustering
    # -------------------------------------------------------------------------
    if 5 in tiers_to_run:
        print(f"[Tier 5] Clustering ({args.cluster_identity*100:.0f}% identity threshold)...")
        t5_output = str(output_dir / "tier5_clustered.fasta")
        start = datetime.now()
        filter_by_clustering(current_input, t5_output, args.cluster_identity)
        log_tier(log, 5, current_input, t5_output, start, datetime.now())
        current_input = t5_output

    # -------------------------------------------------------------------------
    # Tier 6: Molecular Docking
    # -------------------------------------------------------------------------
    if 6 in tiers_to_run:
        if args.structures_dir is None:
            print("[Tier 6] Skipping: --structures_dir not provided (need ESMFold PDB structures)")
        else:
            print(f"[Tier 6] Molecular docking (NADH + DHNA)...")
            t6_output = str(output_dir / "tier6_docking.fasta")
            t6_scores = str(output_dir / "tier6_docking_scores.csv")
            t6_work = str(output_dir / "tier6_docking_work")
            start = datetime.now()
            filter_by_docking(current_input, t6_output, args.structures_dir,
                              args.min_nadh_energy, args.min_dhna_energy,
                              t6_scores, t6_work)
            log_tier(log, 6, current_input, t6_output, start, datetime.now())
            current_input = t6_output

    # -------------------------------------------------------------------------
    # Tier 7: DFT
    # -------------------------------------------------------------------------
    if 7 in tiers_to_run:
        if args.structures_dir is None or args.docking_dir is None:
            print("[Tier 7] Skipping: --structures_dir and --docking_dir required for DFT")
        else:
            print(f"[Tier 7] DFT calculations (HOMO-LUMO gap >= {args.min_gap} eV)...")
            t7_output = str(output_dir / "tier7_dft.fasta")
            t7_scores = str(output_dir / "tier7_dft_scores.csv")
            t7_work = str(output_dir / "tier7_dft_work")
            start = datetime.now()
            filter_by_dft(current_input, t7_output, args.structures_dir,
                          args.docking_dir, args.min_gap, t7_scores, t7_work)
            log_tier(log, 7, current_input, t7_output, start, datetime.now())
            current_input = t7_output

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    log_file = output_dir / "pipeline_log.json"
    with open(log_file, "w") as f:
        json.dump(log, f, indent=2)

    print(f"\n{'='*60}")
    print(f"PIPELINE COMPLETE")
    print(f"{'='*60}")
    print(f"{'Tier':<8} {'Input':>8} {'Output':>8} {'Retained':>10} {'Duration':>12}")
    print(f"{'-'*50}")
    for entry in log:
        retention = 100 * entry["output_count"] / entry["input_count"] if entry["input_count"] > 0 else 0
        print(f"{'Tier '+str(entry['tier']):<8} {entry['input_count']:>8} {entry['output_count']:>8} {retention:>9.1f}% {entry['duration_seconds']:>10.1f}s")
    print(f"{'-'*50}")
    print(f"Final output: {current_input} ({count_seqs(current_input)} sequences)")
    print(f"Log saved to: {log_file}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
