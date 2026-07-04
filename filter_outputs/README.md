# NDH-2 Filtration Pipeline — Organized by Tier

This folder contains scripts and sbatch files for the 7-tier computational filtration
pipeline applied to GenSLM-generated NDH-2 sequences. Each subdirectory corresponds
to a single tier.

## Tier Overview

| Tier | Name | Criterion | Time Scale | Scripts |
|------|------|-----------|------------|---------|
| 1 | Length Filter | 350–800 aa | Seconds | `tier1_length.py`, `filter_generated.py` |
| 2 | BLAST Homology | 30–90 % identity to reference | Minutes | `tier2_blast.py` |
| 3 | ESMFold pLDDT | Mean pLDDT > 70 | Hours (GPU) | `tier3_esmfold.py`, `tier3_esmfold_cli.py` |
| 4 | Membrane Topology | ≥4 TM helices, C-term bias | Hours | `tier4_topology.py`, `tier4_topology_simple.py` |
| 5 | CD-HIT Clustering | 70 % identity threshold | Minutes | `tier5_clustering.py` |
| 6 | Molecular Docking | NADH + DHNA < −7 kcal/mol | Days | `tier6_docking.py` |
| 7 | DFT Validation | HOMO-LUMO gap > 2 eV | Weeks | `tier7_dft.py` |

## Orchestrator

`run_pipeline.py` (top-level in this folder) runs all tiers sequentially.
Usage:
```bash
python filter_outputs/run_pipeline.py \
    --input <input.fasta> \
    --output_dir <results_dir> \
    --tiers 1,2,3,4,5
```

## Folder Structure

```
filter_outputs/
├── README.md                          ← this file
├── run_pipeline.py                    ← orchestrator (all tiers)
├── tier1_length/
│   ├── tier1_length.py                ← core filter (Biopython length check)
│   ├── filter_generated.py            ← codon-level pre-filter (stops, dedup, length)
│   ├── run_pipeline_tiers1_2.sbatch
│   └── run_pipeline_tiers1_2_600max.sbatch
├── tier2_blast/
│   ├── tier2_blast.py                 ← BLAST similarity filter
│   └── run_pipeline_tier2_blast.sbatch
├── tier3_esmfold/
│   ├── tier3_esmfold.py               ← ESMFold pLDDT filter (Python API)
│   ├── tier3_esmfold_cli.py           ← ESMFold via CLI
│   ├── run_pipeline_tier3_esmfold.sbatch
│   ├── run_pipeline_tier3_esmfold_v100.sbatch
│   ├── run_pipeline_tier3_esmfold_cpu.sbatch
│   ├── run_pipeline_tier3_esmfold_cli.sbatch
│   └── run_pipeline_tier3_600max.sbatch
├── tier4_topology/
│   ├── tier4_topology.py              ← DeepTMHMM-based filter
│   ├── tier4_topology_simple.py       ← hydropathy-based fallback
│   ├── tier4_membrane_check.py        ← additional membrane analysis
│   ├── run_pipeline_tier4_topology.sbatch
│   ├── run_pipeline_tier4_topology_simple.sbatch
│   └── run_pipeline_tier4_topology_soluble.sbatch
├── tier5_clustering/
│   └── tier5_clustering.py            ← CD-HIT diversity selection
├── tier6_docking/
│   └── tier6_docking.py               ← AutoDock Vina / GNINA docking
└── tier7_dft/
    └── tier7_dft.py                   ← Gaussian/ORCA DFT calculations
```

## Notes

- **Tiers 1–5** can be run on standard HPC nodes (Tier 3 benefits from GPU).
- **Tier 6** requires AutoDock Vina or GNINA and ESMFold PDB structures from Tier 3.
- **Tier 7** requires Gaussian 16 or ORCA and docking poses from Tier 6.
- The `filter_generated.py` in `tier1_length/` is the codon-level pre-processing step
  (trim stop codons, remove internal stops, AA length 330–600, deduplication) run
  immediately after generation and before formal FASTA-based pipeline tiers.
