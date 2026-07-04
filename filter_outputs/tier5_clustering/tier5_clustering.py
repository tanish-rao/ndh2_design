"""
Tier 5 Filter: Clustering for Diversity
Clusters sequences using CD-HIT and selects representative sequences.
- Removes redundant sequences within the passed set
- Ensures diverse coverage of sequence space
Runs in minutes on ~500 sequences.
"""
import argparse
import subprocess
from pathlib import Path
from Bio import SeqIO


def run_cdhit(input_fasta, output_fasta, identity_threshold=0.7, threads=8):
    """Run CD-HIT to cluster sequences and select representatives."""
    import os
    import shutil
    
    # Try to find cd-hit in conda environment
    cdhit_path = shutil.which("cd-hit")
    if not cdhit_path:
        # Try conda environment path
        conda_env = os.environ.get("CONDA_PREFIX", "")
        if conda_env:
            cdhit_path = os.path.join(conda_env, "bin", "cd-hit")
            if not os.path.exists(cdhit_path):
                cdhit_path = "cd-hit"
        else:
            cdhit_path = "cd-hit"
    
    # Choose word length (-n) based on identity threshold per CD-HIT rules
    # 0.7-1.0 -> n=5, 0.6-0.7 -> n=4, 0.5-0.6 -> n=3, 0.4-0.5 -> n=2
    if identity_threshold >= 0.7:
        n_word = "5"
    elif identity_threshold >= 0.6:
        n_word = "4"
    elif identity_threshold >= 0.5:
        n_word = "3"
    else:
        n_word = "2"

    cmd = [
        cdhit_path,
        "-i", str(input_fasta),
        "-o", str(output_fasta),
        "-c", str(identity_threshold),
        "-n", n_word,
        "-T", str(threads),
        "-M", "16000",
        "-d", "0"
    ]
    print(f"Running CD-HIT (identity threshold: {identity_threshold})...")
    print(f"CD-HIT path: {cdhit_path}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"CD-HIT error: {result.stderr}")
        return False
    return True


def filter_by_clustering(input_fasta, output_fasta, identity_threshold=0.7, threads=8):
    input_count = sum(1 for _ in SeqIO.parse(input_fasta, "fasta"))

    success = run_cdhit(input_fasta, output_fasta, identity_threshold, threads)
    if not success:
        print("CD-HIT failed. Trying MMseqs2 as fallback...")
        success = run_mmseqs2(input_fasta, output_fasta, identity_threshold, threads)
        if not success:
            print("ERROR: Both CD-HIT and MMseqs2 failed. Check installations.")
            return []

    passed = list(SeqIO.parse(output_fasta, "fasta"))
    print(f"Clustering filter ({identity_threshold*100:.0f}% identity threshold):")
    print(f"  Input:  {input_count} sequences")
    print(f"  Output: {len(passed)} representative sequences ({100*len(passed)/input_count:.1f}%)")
    print(f"  Output: {output_fasta}")
    return passed


def run_mmseqs2(input_fasta, output_fasta, identity_threshold=0.7, threads=8):
    """Fallback: MMseqs2 clustering."""
    import tempfile, os
    with tempfile.TemporaryDirectory() as tmpdir:
        db = os.path.join(tmpdir, "seqdb")
        clust = os.path.join(tmpdir, "clust")
        rep = os.path.join(tmpdir, "rep")
        tmp = os.path.join(tmpdir, "tmp")

        subprocess.run(["mmseqs", "createdb", str(input_fasta), db], check=True, capture_output=True)
        subprocess.run(["mmseqs", "cluster", db, clust, tmp,
                        "--min-seq-id", str(identity_threshold),
                        "--threads", str(threads)], check=True, capture_output=True)
        subprocess.run(["mmseqs", "createsubdb", clust, db, rep], check=True, capture_output=True)
        subprocess.run(["mmseqs", "convert2fasta", rep, str(output_fasta)], check=True, capture_output=True)
    return True


def main():
    parser = argparse.ArgumentParser(description="Cluster sequences for diversity selection.")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file (cluster representatives)")
    parser.add_argument("--identity", type=float, default=0.7, help="Clustering identity threshold (default: 0.7)")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads (default: 8)")
    args = parser.parse_args()

    filter_by_clustering(args.input, args.output, args.identity, args.threads)


if __name__ == "__main__":
    main()
