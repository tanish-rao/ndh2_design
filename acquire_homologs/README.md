# Acquire NDH2 Homologs (CDS DNA)

End-to-end pipeline to build a high-confidence NDH2 CDS dataset for Gen-SLM fine-tuning. It fetches CDS DNA via NCBI E-utilities, applies strict off-target exclusion, and confirms candidates by an NDH2 HMM screen.

## Final dataset (Strict+HMM)
- Folder: `acquire_homologs/outputs_strict_hmm/`
  - `ndh2_cds_expanded.fasta` — final CDS DNA sequences
  - `ndh2_cds_expanded_codons.txt` — codon triplets per CDS
  - `ndh2_proteins_expanded.fasta` — translated proteins
  - `ndh2_cds_expanded_metadata.tsv` — nt len, aa len, GC%

## Reproducible pipeline
Requirements:
- Python 3.8+, `requests`
- Slurm cluster with modules: `mafft`, `hmmer`
- Optional: `NCBI_API_KEY` for higher E-utilities rate limits

Steps:
1) Expand and fetch CDS with resume-safe batching
```
sbatch acquire_homologs/run_expand_homologs.sh
```
Outputs to `acquire_homologs/outputs_expanded/`, including `nuccore_ids.txt` and `tmp_cds_batches/` cache.

2) HMMER domain screen to confirm NDH2
```
sbatch acquire_homologs/run_hmmer_screen.sh
```
This builds an NDH2 HMM from a strict seed, prepares AA 330–600 candidates (excluding Complex I terms), runs `hmmsearch`, creates an allowlist of protein_ids, and writes strict+HMM-filtered outputs to `outputs_expanded/`.

3) Strict+HMM snapshot (for training)
The strict+HMM set (AA 330–600, Complex I excluded, domain-confirmed) is copied to `outputs_strict_hmm/` for stable consumption.

## Notes
- The large fetch cache `outputs_expanded/tmp_cds_batches/` can be deleted after producing final outputs (it was used to resume long EFetch jobs).
- HMMER artifacts live in `outputs_expanded/hmmer_screen/` for audit/re-tuning.
- To adjust strictness: edit thresholds in `run_hmmer_screen.sh` or rebuild the allowlist from `hits.domtbl` and rerun `postprocess_expanded_cds.py` with custom parameters.

## Key scripts
- `expand_homologs.py` — Broad ESearch+ELink+EFetch with retry + checkpointing
- `run_expand_homologs.sh` — Slurm wrapper (creates `outputs_expanded/`)
- `postprocess_expanded_cds.py` — Parse, filter, dedup; strict/lenient modes; allowlist support
- `prepare_candidates.py` — Build AA candidates and header map from cached CDS
- `run_hmmer_screen.sh` — MAFFT align strict seeds, build HMM, hmmsearch, allowlist, finalize outputs


## Combined dataset (Strict+HMM ∪ Balanced)
- Folder: `acquire_homologs/combined/`
  - `combined_codon_tokens.txt` — merged codon-token training data (deduplicated by protein)
  - `combined_proteins.fasta` — corresponding proteins
  - `combined_metadata.tsv` — id, source, nt_len, aa_len, GC%, protein k-mer Jaccard to parent (k=3)
  - `combined_train.txt`, `combined_val.txt` — 90/10 split of codon tokens
  - `balanced_hmm.tbl` — HMMER tblout used to filter the prior balanced set

How this was built:
1) Create a protein FASTA from the prior balanced set for HMM screening (auto-generated): `combined/balanced_tmp.faa`.
2) Run HMMER against the existing NDH2 HMM to keep only NDH2-like balanced entries:
```
sbatch acquire_homologs/run_hmm_balanced.sbatch
```
This writes `combined/balanced_hmm.tbl`.
3) Merge with the Strict+HMM codons and deduplicate by exact protein:
```
python -u acquire_homologs/combine_datasets.py \
  ndh2_balanced_for_genslm.txt \
  acquire_homologs/outputs_strict_hmm/ndh2_cds_expanded_codons.txt \
  acquire_homologs/outputs_expanded/hmmer_screen/ndh2.hmm \
  parent_info/parent_enzyme.fasta \
  acquire_homologs/combined \
  acquire_homologs/combined/balanced_hmm.tbl
```

Quick plots for the combined set:
```
python -u acquire_homologs/analyze_combined.py acquire_homologs/combined
```
Outputs to `acquire_homologs/combined/plots/` (AA length, GC%, protein-similarity, length vs GC).

