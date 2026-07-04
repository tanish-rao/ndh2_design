"""
Tier 2 Filter: BLAST Similarity
Filters sequences by BLAST similarity against NDH-2 reference sets.
- Removes sequences with <30% identity (too dissimilar, likely junk)
- Removes sequences with >90% identity (too similar to known sequences, not novel)
Runs in minutes on ~8K sequences.
"""
import argparse
import subprocess
import tempfile
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Blast import NCBIXML


def run_blast(query_fasta, db_fasta, output_xml, evalue=1e-5, num_threads=8):
    """Build BLAST db and run blastp."""
    db_path = str(db_fasta) + "_blastdb"

    print(f"Building BLAST database from {db_fasta}...")
    subprocess.run([
        "makeblastdb", "-in", str(db_fasta), "-dbtype", "prot",
        "-out", db_path, "-parse_seqids"
    ], check=True, capture_output=True)

    print(f"Running BLAST ({query_fasta} vs {db_fasta})...")
    subprocess.run([
        "blastp", "-query", str(query_fasta), "-db", db_path,
        "-out", str(output_xml), "-outfmt", "5",
        "-evalue", str(evalue), "-num_threads", str(num_threads),
        "-max_target_seqs", "5"
    ], check=True)


def parse_blast_results(blast_xml):
    """Parse BLAST XML and return {query_id: best_identity}."""
    results = {}
    with open(blast_xml) as f:
        for record in NCBIXML.parse(f):
            query_id = record.query.split()[0]
            if record.alignments:
                best_hsp = record.alignments[0].hsps[0]
                identity = best_hsp.identities / best_hsp.align_length * 100
                results[query_id] = identity
            else:
                results[query_id] = 0.0
    return results


def filter_by_blast(input_fasta, output_fasta, ref_fasta, min_identity=30.0, max_identity=90.0, num_threads=8):
    sequences = list(SeqIO.parse(input_fasta, "fasta"))

    with tempfile.TemporaryDirectory() as tmpdir:
        blast_xml = os.path.join(tmpdir, "blast_results.xml")
        run_blast(input_fasta, ref_fasta, blast_xml, num_threads=num_threads)
        identity_map = parse_blast_results(blast_xml)

    passed, failed_low, failed_high, no_hit = [], [], [], []
    for record in sequences:
        seq_id = record.id
        identity = identity_map.get(seq_id, 0.0)
        if identity == 0.0:
            no_hit.append(record)
        elif identity < min_identity:
            failed_low.append(record)
        elif identity > max_identity:
            failed_high.append(record)
        else:
            passed.append(record)

    SeqIO.write(passed, output_fasta, "fasta")

    total = len(sequences)
    print(f"BLAST similarity filter ({min_identity}-{max_identity}% identity):")
    print(f"  Input:          {total} sequences")
    print(f"  Passed:         {len(passed)} ({100*len(passed)/total:.1f}%)")
    print(f"  Failed (low):   {len(failed_low)} (<{min_identity}% identity)")
    print(f"  Failed (high):  {len(failed_high)} (>{max_identity}% identity, too similar to known)")
    print(f"  No BLAST hit:   {len(no_hit)} (removed)")
    print(f"  Output: {output_fasta}")
    return passed


def main():
    parser = argparse.ArgumentParser(description="Filter sequences by BLAST similarity.")
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    parser.add_argument("--ref", required=True, help="Reference FASTA file to BLAST against")
    parser.add_argument("--min_identity", type=float, default=30.0, help="Min %% identity (default: 30)")
    parser.add_argument("--max_identity", type=float, default=90.0, help="Max %% identity (default: 90)")
    parser.add_argument("--threads", type=int, default=8, help="Number of BLAST threads (default: 8)")
    args = parser.parse_args()

    filter_by_blast(args.input, args.output, args.ref, args.min_identity, args.max_identity, args.threads)


if __name__ == "__main__":
    main()
