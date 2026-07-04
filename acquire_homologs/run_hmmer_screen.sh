#!/bin/bash
#SBATCH --job-name=ndh2_hmmer
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=acquire_homologs/logs/hmmer_%j.out
#SBATCH --error=acquire_homologs/logs/hmmer_%j.err

set -euo pipefail

OUTDIR=acquire_homologs/outputs_expanded
BATCH_DIR=${OUTDIR}/tmp_cds_batches
WORKDIR=${OUTDIR}/hmmer_screen
mkdir -p "$WORKDIR" acquire_homologs/logs

module purge || true
module load mafft/7.505-gcc-13.2.0-nklkvtc hmmer/3.3.2-gcc-13.2.0-yvuiqw4

# 1) Build seed set from strict-filtered proteins (already written)
SEED_FASTA=${OUTDIR}/ndh2_proteins_expanded.fasta
if [ ! -s "$SEED_FASTA" ]; then
  echo "ERROR: seed fasta not found: $SEED_FASTA" >&2
  exit 1
fi
# Downsample to max 1000 sequences for faster alignment
python - "$SEED_FASTA" > "$WORKDIR/seed_1000.faa" << 'PY'
import sys, random
random.seed(0)
seqs=[]; name=None; parts=[]
for line in open(sys.argv[1]):
    line=line.strip()
    if line.startswith('>'):
        if name is not None:
            seqs.append((name, ''.join(parts)))
        name=line[1:]; parts=[]
    else:
        parts.append(line)
if name is not None:
    seqs.append((name, ''.join(parts)))
random.shuffle(seqs)
seqs=seqs[:1000]
for n,s in seqs:
    print(f'>{n}')
    for i in range(0,len(s),60):
        print(s[i:i+60])
PY

# 2) MSA with MAFFT
mafft --quiet --thread ${SLURM_CPUS_PER_TASK:-8} "$WORKDIR/seed_1000.faa" > "$WORKDIR/seed_1000.aln.faa"

# 3) Build HMM
hmmbuild "$WORKDIR/ndh2.hmm" "$WORKDIR/seed_1000.aln.faa" > "$WORKDIR/hmmbuild.log"

# 4) Prepare candidate DB (AA 330-600, exclude Complex I terms)
python acquire_homologs/prepare_candidates.py "$BATCH_DIR" "$WORKDIR" 330 600

# 5) hmmsearch against candidates
hmmsearch --cpu ${SLURM_CPUS_PER_TASK:-8} --tblout "$WORKDIR/hits.tbl" --domtblout "$WORKDIR/hits.domtbl" \
  "$WORKDIR/ndh2.hmm" "$WORKDIR/candidates.faa" > "$WORKDIR/hmmsearch.log"

# 6) Parse hits and produce allowlist of protein IDs
python - "$WORKDIR/hits.tbl" "$WORKDIR/candidates_map.tsv" > "$WORKDIR/allowlist.txt" << 'PY'
import sys,re
from collections import defaultdict

tbl=sys.argv[1]; mp=sys.argv[2]
keep=set()
# parse tblout (skip comments)
with open(tbl) as f:
    for line in f:
        if not line.strip() or line.startswith('#'): continue
        cols=line.split()
        if len(cols)>=1:
            keep.add(cols[0])  # target name
# map target name back to protein_id
pid=set()
for line in open(mp):
    if line.startswith('name\t'): continue
    name, protein_id, batch, hdr = line.rstrip('\n').split('\t',3)
    if name in keep:
        if protein_id:
            pid.add(protein_id)
# print allowlist
for x in sorted(pid):
    print(x)
PY

# 7) Final strict+HMM filter writing outputs
python acquire_homologs/postprocess_expanded_cds.py "$OUTDIR" "$BATCH_DIR" strict 330 600 "$WORKDIR/allowlist.txt"

echo "HMMER screen complete at $(date)"
