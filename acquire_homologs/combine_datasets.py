#!/usr/bin/env python3
import os
import sys
import re
import random
import subprocess
from collections import OrderedDict

# Simple codon table (NCBI 11 bacterial and plant plastid; standard for NDH2)
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}


def translate_dna(dna):
    dna = dna.upper().replace('U', 'T')
    aa = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        aa.append(CODON_TABLE.get(codon, 'X'))
    prot = ''.join(aa)
    # stop at first stop if present
    stop_idx = prot.find('*')
    if stop_idx != -1:
        prot = prot[:stop_idx]
    return prot


def read_codon_tokens(path, source_label):
    entries = []  # list of dicts: {id, source, dna, prot}
    with open(path, 'r') as f:
        idx = 0
        for line in f:
            line = line.strip()
            if not line or line.startswith('>'):
                continue
            dna = line.replace(' ', '').upper().replace('U', 'T')
            dna = re.sub(r'[^ACGT]', '', dna)
            if not dna:
                continue
            prot = translate_dna(dna)
            entries.append({
                'id': f'{source_label}_{idx}',
                'source': source_label,
                'dna': dna,
                'prot': prot,
            })
            idx += 1
    return entries


def write_fasta(entries, path):
    with open(path, 'w') as f:
        for e in entries:
            f.write(f'>{e["id"]}\n')
            # wrap to 80
            s = e['prot']
            for i in range(0, len(s), 80):
                f.write(s[i:i+80] + '\n')


def kmerset(seq, k):
    return set(seq[i:i+k] for i in range(len(seq) - k + 1)) if len(seq) >= k else set()


def jaccard(a, b):
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def read_fasta_protein(path):
    name = None
    seq = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq)
                name = line[1:].strip()
                seq = []
            else:
                seq.append(re.sub(r'\s+', '', line))
        if name is not None:
            yield name, ''.join(seq)


def gc_pct(dna):
    if not dna:
        return float('nan')
    gc = sum(1 for c in dna if c in ('G','C'))
    return 100.0 * gc / len(dna)


def run_hmmsearch(hmm_path, proteins_faa, tblout_path, evalue=1e-10):
    cmd = [
        'hmmsearch',
        '--noali',
        f'-E', str(evalue),
        '--tblout', tblout_path,
        hmm_path,
        proteins_faa,
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def parse_tblout_ids(tblout_path):
    ids = set()
    with open(tblout_path, 'r') as f:
        for line in f:
            if not line or line.startswith('#'):
                continue
            parts = line.strip().split()
            if len(parts) >= 1:
                target = parts[0]
                ids.add(target)
    return ids


def main():
    # Args: [balanced_codon_txt] [strict_codons_txt] [hmm_path] [parent_protein_fasta] [out_dir] [precomputed_tblout(optional)]
    if len(sys.argv) < 6:
        print('Usage: combine_datasets.py [balanced_codon_txt] [strict_codons_txt] [hmm_path] [parent_protein_fasta] [out_dir] [precomputed_tblout(optional)]')
        sys.exit(1)
    balanced_codon_txt = sys.argv[1]
    strict_codons_txt = sys.argv[2]
    hmm_path = sys.argv[3]
    parent_protein_fasta = sys.argv[4]
    out_dir = sys.argv[5]
    pre_tbl = sys.argv[6] if len(sys.argv) > 6 else ''

    os.makedirs(out_dir, exist_ok=True)

    # Load parent protein for k-mer similarity
    parent_seq = None
    for _, s in read_fasta_protein(parent_protein_fasta):
        parent_seq = s
        break
    parent_kset = kmerset(parent_seq, 3) if parent_seq else set()

    # Load entries
    balanced_entries = read_codon_tokens(balanced_codon_txt, 'balanced')
    strict_entries = read_codon_tokens(strict_codons_txt, 'strict_hmm')

    # Filter by AA length
    AA_MIN, AA_MAX = 330, 600
    balanced_entries = [e for e in balanced_entries if AA_MIN <= len(e['prot']) <= AA_MAX]
    strict_entries = [e for e in strict_entries if AA_MIN <= len(e['prot']) <= AA_MAX]

    # HMMER screen on balanced proteins
    balanced_faa = os.path.join(out_dir, 'balanced_tmp.faa')
    write_fasta(balanced_entries, balanced_faa)
    tblout = os.path.join(out_dir, 'balanced_hmm.tbl')
    allow_ids = set()
    if pre_tbl and os.path.exists(pre_tbl):
        allow_ids = parse_tblout_ids(pre_tbl)
    else:
        try:
            run_hmmsearch(hmm_path, balanced_faa, tblout)
            allow_ids = parse_tblout_ids(tblout)
        except FileNotFoundError as e:
            print('[combine] WARNING: hmmsearch not found; to proceed, run HMMER externally and provide tblout via 6th arg, or install HMMER.')
            raise
    balanced_entries = [e for e in balanced_entries if e['id'] in allow_ids]

    # Deduplicate by exact protein across combined sets (keep first occurrence preferring strict first)
    dedup = OrderedDict()
    for e in strict_entries + balanced_entries:
        key = e['prot']
        if key not in dedup:
            dedup[key] = e
    combined = list(dedup.values())

    # Compute metrics and write outputs
    # Combined codon tokens
    combined_codons_path = os.path.join(out_dir, 'combined_codon_tokens.txt')
    with open(combined_codons_path, 'w') as f:
        for e in combined:
            # space-separated codons
            codons = ' '.join(e['dna'][i:i+3] for i in range(0, len(e['dna']), 3))
            f.write(codons + '\n')

    # Proteins FASTA
    combined_prot_faa = os.path.join(out_dir, 'combined_proteins.fasta')
    write_fasta(combined, combined_prot_faa)

    # Metadata TSV
    meta_path = os.path.join(out_dir, 'combined_metadata.tsv')
    with open(meta_path, 'w') as f:
        f.write('id\tsource\tnt_len\taa_len\tgc_pct\tjaccard_k3_to_parent_protein\n')
        for e in combined:
            nt = len(e['dna'])
            aa = len(e['prot'])
            gc = gc_pct(e['dna'])
            pjac = jaccard(kmerset(e['prot'], 3), parent_kset) if parent_kset else float('nan')
            f.write(f"{e['id']}\t{e['source']}\t{nt}\t{aa}\t{gc:.3f}\t{('%.5f' % pjac) if (pjac==pjac) else ''}\n")

    # Train/val split (90/10)
    rng = random.Random(42)
    idxs = list(range(len(combined)))
    rng.shuffle(idxs)
    split = int(0.9 * len(combined))
    train_idx = set(idxs[:split])
    train_path = os.path.join(out_dir, 'combined_train.txt')
    val_path = os.path.join(out_dir, 'combined_val.txt')
    with open(train_path, 'w') as ft, open(val_path, 'w') as fv:
        for i, e in enumerate(combined):
            codons = ' '.join(e['dna'][j:j+3] for j in range(0, len(e['dna']), 3))
            (ft if i in train_idx else fv).write(codons + '\n')

    print(f"[combine] Strict entries: {len(strict_entries)}; Balanced kept after HMM: {len(balanced_entries)}; Combined (dedup): {len(combined)}")
    print(f"[combine] Wrote: {combined_codons_path}, {combined_prot_faa}, {meta_path}, {train_path}, {val_path}")


if __name__ == '__main__':
    main()
