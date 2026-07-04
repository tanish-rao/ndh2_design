#!/usr/bin/env python3
import os
import sys
import re
from datetime import datetime

def ts():
    return datetime.now().isoformat(timespec='seconds')

GENETIC_CODE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

EXCLUDE_TERMS = [
    'nuod','nuoc','nuob','nuoh','nuoi','nuoj','nuok','nuol','nuom','nuon',
    'complex i','nadh-ubiquinone oxidoreductase','nadh:ubiquinone',
    'nadh dehydrogenase subunit','nadh dehydrogenase i',
]

NDH2_TERMS = [
    # Core gene/product names
    'ndh2','ndh-2','ndh ii','ndh-ii','ndhii',
    'ndi1','ndi-1',
    # Explicit NDH-2 product names
    'nadh dehydrogenase ii',
    'nadh dehydrogenase, type ii',
    'type ii nadh dehydrogenase',
    'type 2 nadh dehydrogenase',
    'type ii nadh:quinone oxidoreductase',
    'nadh:quinone oxidoreductase',
    'nadh dehydrogenase (quinone) 2',
    'nadh dehydrogenase (ubiquinone) 2',
    'nadh dehydrogenase 2',
    # Alternative naming
    'alternative nadh dehydrogenase',
    'alt nadh dehydrogenase',
    # Broader but NDH2-associated
    'nadh dehydrogenase (quinone)',
    'nadh dehydrogenase [quinone]',
    'nadh oxidoreductase',
]

def is_excluded(header: str) -> bool:
    hl = header.lower()
    for t in EXCLUDE_TERMS:
        if t in hl:
            return True
    return False

def is_ndh2_like(header: str) -> bool:
    hl = header.lower()
    for t in NDH2_TERMS:
        if t in hl:
            return True
    if 'nadh dehydrogenase' in hl and not is_excluded(header):
        return True
    return False

def is_ndh2_strict(header: str) -> bool:
    """Strict NDH2 call using explicit synonyms or structured fields.
    Accept if:
      - header contains known NDH2 synonyms (NDH-2, NDH2, NADH dehydrogenase II, etc.), OR
      - [gene=] is ndh/ndh2/ndi1, OR
      - [protein=] contains NADH dehydrogenase/oxidoreductase but NOT 'subunit' and not Complex I markers.
    """
    hl = header.lower()
    # Quick reject if excluded terms present
    if is_excluded(header):
        return False
    # 1) Synonym hits
    for t in NDH2_TERMS:
        if t in hl:
            return True
    # 2) Parse structured fields
    gene = None
    protein = None
    mg = re.search(r'gene=([^\]]+)', header, flags=re.IGNORECASE)
    if mg:
        gene = mg.group(1).strip().lower()
    mp = re.search(r'protein=([^\]]+)', header, flags=re.IGNORECASE)
    if mp:
        protein = mp.group(1).strip().lower()
    if gene in { 'ndh', 'ndh2', 'ndi1' }:
        return True
    if protein:
        if ('nadh dehydrogenase' in protein or 'quinone oxidoreductase' in protein or 'nadh oxidoreductase' in protein):
            # Exclude Complex I wording
            if any(x in protein for x in ['subunit', 'complex i', 'dehydrogenase i', 'nuo']):
                return False
            return True
    return False

def translate(seq: str) -> str:
    aa = []
    for i in range(0, len(seq)-2, 3):
        aa.append(GENETIC_CODE.get(seq[i:i+3], 'X'))
    return ''.join(aa)

def valid_cds(seq: str) -> bool:
    if len(seq) < 3 or len(seq) % 3 != 0:
        return False
    if not all(c in 'ATGC' for c in seq):
        return False
    return True

def fasta_stream(path):
    with open(path) as f:
        header = None
        parts = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(parts).upper().replace(' ', '')
                header = line[1:].strip()
                parts = []
            else:
                parts.append(line)
        if header is not None:
            yield header, ''.join(parts).upper().replace(' ', '')

def process_entries(stream_paths, out_dir, prev_merge_path=None, min_aa=330, max_aa=500, strict=False, protein_ids=None):
    os.makedirs(out_dir, exist_ok=True)
    out_fasta = open(os.path.join(out_dir, 'ndh2_cds_expanded.fasta'), 'w')
    out_codons = open(os.path.join(out_dir, 'ndh2_cds_expanded_codons.txt'), 'w')
    out_meta = open(os.path.join(out_dir, 'ndh2_cds_expanded_metadata.tsv'), 'w')
    out_prot = open(os.path.join(out_dir, 'ndh2_proteins_expanded.fasta'), 'w')
    out_meta.write('header\tnt_len\taa_len\tgc_pct\n')

    seen = set()
    kept = 0
    total = 0

    def maybe_write(h, s):
        nonlocal kept
        if s in seen:
            return
        seen.add(s)
        aa = translate(s)
        if aa.endswith('*'):
            aa = aa[:-1]
        out_fasta.write(f'>{h}\n')
        for i in range(0, len(s), 60):
            out_fasta.write(s[i:i+60] + '\n')
        out_codons.write(f'>{h}\n')
        out_codons.write(' '.join(s[i:i+3] for i in range(0, len(s), 3)) + '\n')
        out_prot.write(f'>{h}\n')
        for i in range(0, len(aa), 60):
            out_prot.write(aa[i:i+60] + '\n')
        gc = 100.0 * (s.count('G') + s.count('C')) / max(1, len(s))
        out_meta.write(f'{h}\t{len(s)}\t{len(s)//3}\t{gc:.2f}\n')
        kept += 1

    def consider(h, s):
        # If a protein ID allowlist is provided, ensure header matches one of them
        if protein_ids is not None and len(protein_ids) > 0:
            # Common patterns in fasta_cds_na headers include [protein_id=ACC]
            m = re.search(r'protein_id=([^\] ]+)', h)
            if not m:
                return
            pid = m.group(1)
            if pid not in protein_ids:
                return
        if is_excluded(h):
            return
        if not valid_cds(s):
            return
        aa_len = len(s) // 3
        if aa_len < min_aa or aa_len > max_aa:
            return
        aa = translate(s)
        if aa.count('X') > 5:
            return
        if aa and aa[0] != 'M' and 'M' not in aa[1:3]:
            return
        maybe_write(h, s)

    # Merge previous filtered set first (if provided)
    if prev_merge_path and os.path.exists(prev_merge_path):
        for h, s in fasta_stream(prev_merge_path):
            total += 1
            consider(h, s)
        print(f'[{ts()}] Merged previous set: {kept} kept so far')

    # Process streamed batch files
    for p in stream_paths:
        for h, s in fasta_stream(p):
            total += 1
            # Strict: require explicit NDH2-like header; Lenient: NDH2-like or non-excluded generic
            if strict:
                if is_ndh2_strict(h):
                    consider(h, s)
            else:
                if is_ndh2_like(h) or not is_excluded(h):
                    consider(h, s)

    for fh in (out_fasta, out_codons, out_meta, out_prot):
        fh.close()

    # Summary
    with open(os.path.join(out_dir, 'expansion_summary.txt'), 'w') as sf:
        sf.write(f'Total entries scanned: {total}\n')
        sf.write(f'Total unique CDS sequences: {kept}\n')
        sf.write(f'Outputs in: {out_dir}/\n')
    print(f'[{ts()}] Wrote summary. Total kept: {kept}')


def main():
    # Args: out_dir [batch_dir] [strict|lenient] [min_aa] [max_aa]
    out_dir = sys.argv[1] if len(sys.argv) > 1 else 'acquire_homologs/outputs_expanded'
    # Determine batch_dir: explicit arg or default inside out_dir; if missing, try outputs_expanded
    if len(sys.argv) > 2 and sys.argv[2] not in ('strict', 'lenient'):
        batch_dir = sys.argv[2]
        arg_idx = 3
    else:
        batch_dir = os.path.join(out_dir, 'tmp_cds_batches')
        arg_idx = 2
        if not os.path.isdir(batch_dir):
            fallback = 'acquire_homologs/outputs_expanded/tmp_cds_batches'
            if os.path.isdir(fallback):
                batch_dir = fallback
    strict = False
    if len(sys.argv) > arg_idx:
        mode = sys.argv[arg_idx].lower()
        if mode == 'strict':
            strict = True
        elif mode == 'lenient':
            strict = False
        else:
            # mode omitted; shift indices
            arg_idx -= 1
    min_aa = int(sys.argv[arg_idx+1]) if len(sys.argv) > arg_idx+1 else 330
    max_aa = int(sys.argv[arg_idx+2]) if len(sys.argv) > arg_idx+2 else 500
    protein_ids_file = sys.argv[arg_idx+3] if len(sys.argv) > arg_idx+3 else ''
    protein_ids = None
    if protein_ids_file and os.path.exists(protein_ids_file):
        with open(protein_ids_file) as pf:
            protein_ids = set(x.strip() for x in pf if x.strip())

    if not os.path.isdir(batch_dir):
        print(f'ERROR: batch dir not found: {batch_dir}')
        sys.exit(1)
    # Collect batch files in order
    batch_files = sorted([os.path.join(batch_dir, x) for x in os.listdir(batch_dir) if re.match(r'^batch_\d+\.fa$', x)])
    prev_path = 'acquire_homologs/outputs/homologs_filtered.fasta'
    print(f'[{ts()}] Starting postprocess on {len(batch_files)} batches (strict={strict}, aa {min_aa}-{max_aa}, protein_ids={"yes" if protein_ids else "no"})')
    process_entries(batch_files, out_dir, prev_merge_path=prev_path, min_aa=min_aa, max_aa=max_aa, strict=strict, protein_ids=protein_ids)
    print(f'[{ts()}] Done')

if __name__ == '__main__':
    main()
