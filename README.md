# NDH2 Design: End-to-End Pipeline

This repository contains an end-to-end pipeline for curating NDH-2 homologs, training a sequence language model, generating candidate sequences focused on the catalytic domain, and filtering candidates through a multi-tier screening workflow to identify high-confidence, diverse designs suitable for downstream analysis and fine-tuning.

## Pipeline at a Glance

```mermaid
flowchart TD
  A[Acquire NDH2 homologs<br/>(Strict + HMM + balanced)] --> B[Combine + QC + Deduplicate]
  B --> C[Train Gen-SLM on catalytic domain]
  C --> D[Generate codon sequences]
  D --> E[Pre-clean codons<br/>trim stops, no internal stops,<br/>AA len bounds, dedupe]
  E --> T1[Tier 1: Length filter]
  T1 --> T2[Tier 2: BLAST filter]
  T2 --> L[Append parent linker (AA 423+)]
  L --> T3[Tier 3: ESMFold pLDDT]
  T3 --> P[Compare pre vs post-linker pLDDT<br/>select top-N by pLDDT]
  P --> T5[Tier 5: Clustering (CD-HIT)]
  T5 --> T6[Tier 6: Docking NADH & DHNA]
  T6 --> T6_1[[Step 6.1: Dock FAD (optional)]]
  T6 --> T7[Tier 7: DFT (optional)]
```

- Linker appending uses the parent enzyme AA 423 onward to retain the catalytic domain focus.
- Tier 4 (topology) is usually skipped for linker-appended constructs (TM helices reside in the linker region), but tooling is available if needed.
- Tier 5 supports multiple identity thresholds (e.g., 70% default, 50–60% supported), with automatic CD-HIT word-length selection.
- Tier 6 filters on NADH and DHNA; Step 6.1 adds FAD docking for pocket profiling (not gating by default).

## Repository Structure

- parent_info/ — Parent enzyme FASTA and metadata
- acquire_homologs/ — HMMER/BLAST-based acquisition and re-screening utilities
- train_model/ — Training and sequence generation utilities
- filter_outputs/
  - run_pipeline.py — Orchestrator for tiered filtering
  - tier1_length/, tier2_blast/, tier3_esmfold/, tier4_topology/, tier5_clustering/, tier6_docking/, tier7_dft/
  - pipeline_results/ — Canonical outputs per tier (FASTA, CSV, figures, logs)
- pipeline_results/ — Prior runs and comparative analyses (historical)

Key outputs live under filter_outputs/pipeline_results/:
- tier1_length.fasta, tier2_blast.fasta, tier3_esmfold.fasta
- tier3_plddt_scores.csv, figs/
- tier5_clustered.fasta (+ .clstr)
- tier6_docking.fasta, tier6_docking_scores.csv, tier6_docking_work/
- pipeline_log.json

## How to Run

Below are the typical steps for the generated-sequence filtration path. Paths assume repository root.

1) Pre-clean generated codons and convert to proteins (already scripted in filter_outputs)
- Input combined/cleaned codons written to filter_outputs/pipeline_results/input_generated_cleaned.fasta

2) Tier 1 and Tier 2 (length + BLAST)
- Outputs into filter_outputs/pipeline_results/

3) Append parent linker (AA 423+), then Tier 3 ESMFold pLDDT
- Pre- and post-linker pLDDT can be compared via compare_plddt_pre_post.py

4) Select top-N by pLDDT (post-linker), then Tier 5 clustering
- 70% identity default; 50–60% supported. CD-HIT word length auto-adjusted.

Example (Tier 5 at 70%):

```bash
sbatch --export=ALL,INPUT_FASTA=filter_outputs/pipeline_results/tier3_esmfold.fasta,IDENT=0.7 \
  filter_outputs/tier5_clustering/run_tier5_on_postlinker.sbatch
```

5) Generate structures for Tier 5 reps (if not already present)

- ESMFold structures are expected at:
  - filter_outputs/pipeline_results/tier5_structures/{seq_id}.pdb

6) Docking (Tier 6) with NADH and DHNA; add Step 6.1 FAD docking (scores only)

```bash
python filter_outputs/tier6_docking/tier6_docking.py \
  --input filter_outputs/pipeline_results/tier5_clustered.fasta \
  --structures_dir filter_outputs/pipeline_results/tier5_structures \
  --output filter_outputs/pipeline_results/tier6_docking.fasta \
  --scores_file filter_outputs/pipeline_results/tier6_docking_scores.csv \
  --work_dir filter_outputs/pipeline_results/tier6_docking_work \
  --include_fad
```

- Filtering is based on NADH and DHNA thresholds (default both ≤ -7.0 kcal/mol). FAD is recorded for analysis but not gating by default. You can override thresholds via --min_nadh_energy and --min_dhna_energy.

7) Optional DFT (Tier 7)
- Tier 7 scripts provide optional quantum refinement on top candidates.

## Figures and Reporting

- Figures are saved to filter_outputs/pipeline_results/figs
- A per-run summary is available at filter_outputs/pipeline_results/pipeline_log.json

## Dependencies (typical)

- Python 3.9+
- PyTorch, Biopython, esm
- CD-HIT, MMseqs2
- Open Babel (obabel)
- AutoDock Vina or GNINA
- DeepTMHMM (if running Tier 4)

Consider managing via conda and a cluster scheduler (Slurm) for GPU/CPU partitioning.

## Notes

- Identity thresholds in clustering refer to sequence identity used to collapse near-duplicates. Lower identity → fewer representatives (more merging), higher identity → more representatives.
- The linker-appended constructs carry transmembrane helices; topology filtering (Tier 4) may be skipped accordingly.
- The new Step 6.1 FAD docking is intended to profile the FAD-binding pocket and complement NADH/DHNA docking.

## Citation

- Lin et al., ESMFold: End-to-end single-sequence structure prediction.
- Li & Godzik, CD-HIT: a fast program for clustering.
- Steinegger & Söding, MMseqs2.
- Trott & Olson, AutoDock Vina.
- Hall, Open Babel.

## Pipeline (ASCII fallback)

```
[Acquire NDH2 homologs] -> [Combine/QC/Dedup] -> [Train Gen-SLM (catalytic)]
    -> [Generate codons] -> [Pre-clean codons]
    -> [Tier 1: Length] -> [Tier 2: BLAST]
    -> [Append parent linker (AA 423+)] -> [Tier 3: ESMFold pLDDT]
    -> [Compare pre/post pLDDT + Select top-N]
    -> [Tier 5: Clustering (CD-HIT)]
    -> [Tier 6: Dock NADH & DHNA]
    -> [Step 6.1: Dock FAD (optional, scoring)]
    -> [Tier 7: DFT (optional)]
```

### Step-by-step (quick reference)

1. Acquire homologs (strict+HMM+balanced), combine, QC, deduplicate.
2. Train Gen-SLM on catalytic domain-focused set; generate codon sequences.
3. Pre-clean codons: trim trailing stops, remove internal stops, enforce AA length, deduplicate.
4. Tier 1: Length filter; Tier 2: BLAST filter.
5. Append parent linker domain (AA 423+), preserving catalytic domain focus.
6. Tier 3: ESMFold pLDDT on linker-appended sequences; optionally compare pre/post.
7. Select top-N by pLDDT (post-linker) for clustering.
8. Tier 5: Cluster (CD-HIT) at chosen identity (default 70%; supports 50–60% with auto word-length).
9. Generate structures for Tier 5 representatives (ESMFold PDBs) if not already present.
10. Tier 6: Dock NADH and DHNA (filtering criteria; default ≤ -7.0 kcal/mol each).
11. Step 6.1 (optional): Dock FAD to score the cofactor pocket (not gating by default).
12. Tier 7 (optional): DFT refinement on top candidates.

