#!/usr/bin/env python3
import os
import json
import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO


def count_fasta(path: Path) -> int:
    if not path or not path.exists():
        return 0
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def lengths_from_fasta(path: Path) -> List[int]:
    if not path or not path.exists():
        return []
    return [len(rec.seq) for rec in SeqIO.parse(str(path), "fasta")]


def try_load_json(path: Path):
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return None


def plot_counts_by_tier(results_dir: Path, tiers_present: Dict[int, Path], log: List[Dict]):
    labels = []
    counts = []

    # If log exists, use that ordering
    if log:
        for entry in log:
            labels.append(f"T{entry['tier']}")
            counts.append(entry["output_count"])  # count after tier
    else:
        # Fallback: compute from discovered files
        for tier in sorted(tiers_present.keys()):
            labels.append(f"T{tier}")
            counts.append(count_fasta(tiers_present[tier]))

    if not counts:
        return

    fig = plt.figure(figsize=(6, 4))
    plt.bar(labels, counts, color="#4e79a7")
    plt.ylabel("Sequences")
    plt.title("Retention by Tier")
    plt.tight_layout()
    fig.savefig(results_dir / "figs" / "counts_by_tier.png", dpi=150)
    plt.close(fig)


def plot_length_hists(results_dir: Path, tiers_present: Dict[int, Path]):
    for tier, fasta in sorted(tiers_present.items()):
        lens = lengths_from_fasta(fasta)
        if not lens:
            continue
        fig = plt.figure(figsize=(6, 4))
        plt.hist(lens, bins=40, color="#59a14f")
        plt.xlabel("AA length")
        plt.ylabel("Count")
        plt.title(f"Tier {tier} length distribution (n={len(lens)})")
        plt.tight_layout()
        fig.savefig(results_dir / "figs" / f"tier{tier}_length_hist.png", dpi=150)
        plt.close(fig)


def plot_plddt_hist(results_dir: Path):
    csv_path = results_dir / "tier3_plddt_scores.csv"
    if not csv_path.exists():
        return
    try:
        import csv
        vals = []
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                v = row.get("plddt")
                if v is not None and str(v).strip() != "":
                    try:
                        vals.append(float(v))
                    except ValueError:
                        pass
        if not vals:
            return
        fig = plt.figure(figsize=(6, 4))
        plt.hist(vals, bins=40, color="#f28e2b")
        plt.xlabel("pLDDT")
        plt.ylabel("Count")
        plt.title("Tier 3 pLDDT distribution")
        plt.tight_layout()
        fig.savefig(results_dir / "figs" / "tier3_plddt_hist.png", dpi=150)
        plt.close(fig)
    except Exception as e:
        print(f"[warn] Failed to plot pLDDT histogram: {e}")


def plot_topology_hist(results_dir: Path):
    # Prefer predictions CSV if present (from simple topology script)
    csv_path = results_dir / "tier4_topology_predictions.csv"
    tm_counts: List[int] = []

    if csv_path.exists():
        try:
            import csv
            with open(csv_path) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    v = row.get("tm_count")
                    if v is not None and str(v).strip() != "":
                        try:
                            tm_counts.append(int(float(v)))
                        except ValueError:
                            pass
        except Exception as e:
            print(f"[warn] Failed to parse topology predictions CSV: {e}")

    # If not found, try parsing GFF3 directory produced by DeepTMHMM fallback
    if not tm_counts:
        gff_dir = results_dir / "tier4_deeptmhmm"
        if gff_dir.exists():
            for gff in gff_dir.glob("*.gff3"):
                try:
                    with open(gff) as f:
                        tm_for_record = 0
                        for line in f:
                            if line.startswith("##"):
                                continue
                            if line.strip() == "//":
                                if tm_for_record > 0:
                                    tm_counts.append(tm_for_record)
                                tm_for_record = 0
                                continue
                            parts = line.split("\t")
                            if len(parts) > 2 and parts[2] == "TMhelix":
                                tm_for_record += 1
                except Exception:
                    pass

    if not tm_counts:
        return

    fig = plt.figure(figsize=(6, 4))
    plt.hist(tm_counts, bins=range(0, max(tm_counts) + 2), align="left", color="#e15759")
    plt.xlabel("Predicted TM helices")
    plt.ylabel("Count")
    plt.title("Tier 4 TM helix count distribution")
    plt.tight_layout()
    fig.savefig(results_dir / "figs" / "tier4_tm_count_hist.png", dpi=150)
    plt.close(fig)


def discover_tier_fastas(results_dir: Path) -> Dict[int, Path]:
    mapping = {}
    # The orchestrator writes canonical names; also support any tierX_*.fasta
    for tier in range(1, 8):
        candidates = list(results_dir.glob(f"tier{tier}_*.fasta"))
        if candidates:
            # Prefer the canonical names if present
            canon = results_dir / {
                1: "tier1_length.fasta",
                2: "tier2_blast.fasta",
                3: "tier3_esmfold.fasta",
                4: "tier4_topology.fasta",
                5: "tier5_clustered.fasta",
                6: "tier6_docking.fasta",
                7: "tier7_dft.fasta",
            }[tier]
            mapping[tier] = canon if canon.exists() else candidates[0]
    return mapping


def main():
    ap = argparse.ArgumentParser(description="Generate plots for NDH-2 pipeline outputs")
    ap.add_argument("--results_dir", default="pipeline_results", help="Directory with pipeline outputs")
    args = ap.parse_args()

    results_dir = Path(args.results_dir)
    figs_dir = results_dir / "figs"
    figs_dir.mkdir(parents=True, exist_ok=True)

    log = try_load_json(results_dir / "pipeline_log.json") or []
    tiers_present = discover_tier_fastas(results_dir)

    plot_counts_by_tier(results_dir, tiers_present, log)
    plot_length_hists(results_dir, tiers_present)
    plot_plddt_hist(results_dir)
    plot_topology_hist(results_dir)

    print(f"Saved figures to: {figs_dir}")


if __name__ == "__main__":
    main()
