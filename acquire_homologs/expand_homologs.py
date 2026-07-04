#!/usr/bin/env python3
"""
Expand NDH2 homolog dataset via:
1. NCBI BLAST of parent protein against nr (bacteria), up to 5000 hits
2. Broader ESearch queries (multiple term variants)
3. Combine protein IDs, link to nuccore, fetch CDS DNA
4. Translate, filter by length/identity, exclude Complex I (nuo) hits
5. Deduplicate
"""
import os
import sys
import time
import re
import json
from typing import List, Tuple, Set
from datetime import datetime
from collections import Counter

try:
    import requests
    from requests.exceptions import ChunkedEncodingError, ConnectionError, ReadTimeout
except Exception:
    print('ERROR: requires requests. pip install requests')
    sys.exit(1)

NCBI_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
BLAST_URL = 'https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi'
UA = {'User-Agent': 'ndh2-homolog-expander/1.0'}

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
    'nuod', 'nuoc', 'nuob', 'nuoh', 'nuoi', 'nuoj', 'nuok', 'nuol', 'nuom', 'nuon',
    'complex i', 'nadh-ubiquinone oxidoreductase', 'nadh:ubiquinone',
    'nadh dehydrogenase subunit', 'nadh dehydrogenase i',
]

def ts():
    return datetime.now().isoformat(timespec='seconds')

def _fetch(url, params, timeout=120, retries=6, as_json=False):
    delay = 2.0
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, headers=UA, timeout=timeout)
            r.raise_for_status()
            if as_json:
                # Clean control characters before parsing
                txt = r.text
                txt = re.sub(r'[\x00-\x08\x0b\x0c\x0e-\x1f]', '', txt)
                return json.loads(txt)
            return r.text
        except (ChunkedEncodingError, ConnectionError, ReadTimeout, json.JSONDecodeError) as e:
            if attempt == retries - 1:
                raise
            print(f"  Retry {attempt+1}/{retries} after {type(e).__name__}")
            time.sleep(delay)
            delay *= 2

def _post(url, data, timeout=120, retries=4):
    delay = 2.0
    for attempt in range(retries):
        try:
            r = requests.post(url, data=data, headers=UA, timeout=timeout)
            r.raise_for_status()
            return r.text
        except (ChunkedEncodingError, ConnectionError, ReadTimeout) as e:
            if attempt == retries - 1:
                raise
            time.sleep(delay)
            delay *= 2

# ─────────────────────────────────────────────────────────────────────
# 1. BLAST parent protein against nr (bacteria)
# ─────────────────────────────────────────────────────────────────────
def submit_blast(query_seq: str, max_hits: int = 5000) -> str:
    params = {
        'CMD': 'Put',
        'PROGRAM': 'blastp',
        'DATABASE': 'nr',
        'QUERY': query_seq,
        'ENTREZ_QUERY': 'bacteria[Organism]',
        'EXPECT': '1e-20',
        'MAX_NUM_SEQ': str(max_hits),
        'FORMAT_TYPE': 'JSON2',
        'HITLIST_SIZE': str(max_hits),
    }
    print(f"[{ts()}] Submitting BLAST job (max {max_hits} hits)...")
    text = _post(BLAST_URL, data=params, timeout=120)
    m = re.search(r'RID\s*=\s*(\S+)', text)
    if not m:
        print('ERROR: Could not parse RID from BLAST submission')
        print(text[:500])
        sys.exit(1)
    rid = m.group(1)
    print(f"[{ts()}] BLAST RID: {rid}")
    return rid

def poll_blast(rid: str, max_wait: int = 3600) -> str:
    start = time.time()
    while time.time() - start < max_wait:
        time.sleep(30)
        params = {
            'CMD': 'Get',
            'RID': rid,
            'FORMAT_TYPE': 'Text',
            'DESCRIPTIONS': '10000',
            'ALIGNMENTS': '0',
        }
        text = _fetch(BLAST_URL, params, timeout=180)
        if 'Status=WAITING' in text:
            elapsed = int(time.time() - start)
            print(f"  [{ts()}] BLAST still running... ({elapsed}s elapsed)")
            continue
        if 'Status=FAILED' in text:
            print('ERROR: BLAST job failed')
            sys.exit(1)
        if 'Status=UNKNOWN' in text:
            print('ERROR: BLAST RID expired or unknown')
            sys.exit(1)
        # Results ready
        return text
    print('ERROR: BLAST timed out')
    sys.exit(1)

def parse_blast_accessions(text: str) -> List[str]:
    """Parse accessions from BLAST Text output (description lines)."""
    accs = set()
    # Match accession patterns like WP_012345678.1, ABC12345.1, etc.
    for m in re.finditer(r'\b([A-Z]{2,3}_?\d{5,12}\.\d+)\b', text):
        accs.add(m.group(1))
    # Also try gi-style: gi|12345|ref|WP_xxx.1|
    for m in re.finditer(r'ref\|([A-Z]{2,3}_?\d+\.\d+)\|', text):
        accs.add(m.group(1))
    for m in re.finditer(r'gb\|([A-Z]{2,3}\d+\.\d+)\|', text):
        accs.add(m.group(1))
    return list(accs)

# ─────────────────────────────────────────────────────────────────────
# 2. Broader ESearch queries
# ─────────────────────────────────────────────────────────────────────
def esearch_protein_broad(api_key: str = '') -> Set[str]:
    queries = [
        '((ndh2[Gene Name]) OR (ndh-2[All Fields]) OR ("NADH dehydrogenase (quinone)"[Title]) OR ("type II NADH dehydrogenase"[Title]) OR ("NDH-2"[Title]) OR ("NADH:quinone oxidoreductase"[Title])) AND bacteria[Organism] NOT hypothetical NOT partial',
        '("type II NADH dehydrogenase"[All Fields]) AND bacteria[Organism]',
        '("NADH dehydrogenase 2"[All Fields]) AND bacteria[Organism]',
        '("NDH-2"[All Fields]) AND bacteria[Organism] NOT "complex I"',
        '("ndh2"[Gene Name]) AND bacteria[Organism]',
        '("ndi1"[Gene Name]) AND bacteria[Organism]',
    ]
    all_ids = set()
    for q in queries:
        params = {
            'db': 'protein',
            'retmode': 'json',
            'retmax': '5000',
            'term': q,
        }
        if api_key:
            params['api_key'] = api_key
        data = _fetch(f'{NCBI_BASE}/esearch.fcgi', params, timeout=60, as_json=True)
        ids = data.get('esearchresult', {}).get('idlist', [])
        all_ids.update(ids)
        print(f"  [{ts()}] ESearch '{q[:60]}...' -> {len(ids)} IDs (cumulative {len(all_ids)})")
        time.sleep(0.5)
    return all_ids

# ─────────────────────────────────────────────────────────────────────
# 3. Convert accessions to protein GIs/IDs
# ─────────────────────────────────────────────────────────────────────
def accessions_to_ids(accs: List[str], api_key: str = '') -> List[str]:
    ids = set()
    B = 200
    for i in range(0, len(accs), B):
        batch = accs[i:i+B]
        params = {
            'db': 'protein',
            'retmode': 'json',
            'retmax': str(len(batch)),
            'term': ' OR '.join(f'{a}[Accession]' for a in batch),
        }
        if api_key:
            params['api_key'] = api_key
        data = _fetch(f'{NCBI_BASE}/esearch.fcgi', params, timeout=60, as_json=True)
        batch_ids = data.get('esearchresult', {}).get('idlist', [])
        ids.update(batch_ids)
        time.sleep(0.5)
    return list(ids)

# ─────────────────────────────────────────────────────────────────────
# 4. Link protein -> nuccore
# ─────────────────────────────────────────────────────────────────────
def elink_to_nuccore(prot_ids: List[str], api_key: str = '') -> List[str]:
    """Link protein IDs to nuccore using XML (more reliable than JSON)."""
    nuccore = set()
    B = 50  # small batches to avoid malformed responses
    total = len(prot_ids)
    n_batches = (total + B - 1) // B
    for i in range(0, total, B):
        batch = prot_ids[i:i+B]
        params = {
            'dbfrom': 'protein',
            'db': 'nuccore',
            'id': ','.join(batch),
            'retmode': 'xml',
        }
        if api_key:
            params['api_key'] = api_key
        try:
            xml_text = _fetch(f'{NCBI_BASE}/elink.fcgi', params, timeout=90)
            # Parse Link IDs from XML with regex (robust)
            for m in re.finditer(r'<Link>\s*<Id>(\d+)</Id>', xml_text):
                nuccore.add(m.group(1))
        except Exception as e:
            print(f"  WARNING: elink batch {i//B+1} failed: {e}")
        batch_num = i // B + 1
        if batch_num % 20 == 0 or batch_num == n_batches:
            print(f"  [{ts()}] elink batch {batch_num}/{n_batches}, nuccore so far: {len(nuccore)}")
        time.sleep(0.4)
    return list(nuccore)

# ─────────────────────────────────────────────────────────────────────
# 5. Fetch CDS DNA
# ─────────────────────────────────────────────────────────────────────
def efetch_cds_na(nuccore_ids: List[str], out_dir: str, api_key: str = '') -> str:
    """Fetch CDS DNA with on-disk checkpointing per batch to allow resume."""
    fasta = ''
    B = 20
    total = len(nuccore_ids)
    n_batches = (total + B - 1) // B
    tmp_dir = os.path.join(out_dir, 'tmp_cds_batches')
    os.makedirs(tmp_dir, exist_ok=True)
    for i in range(0, total, B):
        batch_num = i // B + 1
        batch = nuccore_ids[i:i+B]
        batch_path = os.path.join(tmp_dir, f'batch_{batch_num:04d}.fa')
        # Skip if already fetched and non-empty
        if os.path.exists(batch_path) and os.path.getsize(batch_path) > 0:
            if batch_num % 25 == 0 or batch_num == n_batches:
                print(f"  [{ts()}] efetch batch {batch_num}/{n_batches} (cached)")
            continue
        params = {
            'db': 'nuccore',
            'rettype': 'fasta_cds_na',
            'retmode': 'text',
            'id': ','.join(batch),
        }
        if api_key:
            params['api_key'] = api_key
        text = _fetch(f'{NCBI_BASE}/efetch.fcgi', params, timeout=180)
        with open(batch_path, 'w') as bf:
            bf.write(text)
        if batch_num % 20 == 0 or batch_num == n_batches:
            print(f"  [{ts()}] efetch batch {batch_num}/{n_batches}")
        time.sleep(0.4)
    # Concatenate all batch files
    parts = []
    for b in range(1, n_batches + 1):
        p = os.path.join(tmp_dir, f'batch_{b:04d}.fa')
        if os.path.exists(p):
            with open(p) as f:
                parts.append(f.read())
    return ''.join(parts)

# ─────────────────────────────────────────────────────────────────────
# 6. Parse, translate, filter
# ─────────────────────────────────────────────────────────────────────
def parse_fasta(text: str) -> List[Tuple[str, str]]:
    entries = []
    header = None
    seq_lines = []
    for line in text.splitlines():
        if not line:
            continue
        if line.startswith('>'):
            if header is not None:
                entries.append((header, ''.join(seq_lines).upper().replace(' ', '')))
            header = line[1:].strip()
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header is not None:
        entries.append((header, ''.join(seq_lines).upper().replace(' ', '')))
    return entries

def translate(dna: str) -> str:
    aa = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        aa.append(GENETIC_CODE.get(codon, 'X'))
    return ''.join(aa)

def is_valid_cds(seq: str) -> bool:
    if len(seq) < 3:
        return False
    if len(seq) % 3 != 0:
        return False
    if not all(c in 'ATGC' for c in seq):
        return False
    return True

def is_excluded(header: str) -> bool:
    hl = header.lower()
    for term in EXCLUDE_TERMS:
        if term in hl:
            return True
    return False

def is_ndh2_like(header: str) -> bool:
    hl = header.lower()
    ndh2_terms = [
        'ndh2', 'ndh-2', 'ndi1', 'ndi-1',
        'nadh dehydrogenase (quinone)',
        'type ii nadh dehydrogenase',
        'type 2 nadh dehydrogenase',
        'nadh:quinone oxidoreductase',
        'nadh dehydrogenase 2',
        'nadh oxidoreductase',
        'nadh dehydrogenase [quinone]',
    ]
    for t in ndh2_terms:
        if t in hl:
            return True
    # Also accept generic "NADH dehydrogenase" if not excluded
    if 'nadh dehydrogenase' in hl and not is_excluded(header):
        return True
    return False

def filter_sequences(entries: List[Tuple[str, str]], min_aa=330, max_aa=500) -> List[Tuple[str, str]]:
    passed = []
    for h, s in entries:
        if is_excluded(h):
            continue
        if not is_valid_cds(s):
            continue
        aa_len = len(s) // 3
        if aa_len < min_aa or aa_len > max_aa:
            continue
        protein = translate(s)
        if protein.count('X') > 5:
            continue
        # Must start with M (or near-start)
        if protein[0] != 'M' and protein[1:3].count('M') == 0:
            continue
        passed.append((h, s))
    return passed

def dedup_sequences(entries: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
    seen = set()
    out = []
    for h, s in entries:
        if s in seen:
            continue
        seen.add(s)
        out.append((h, s))
    return out

# ─────────────────────────────────────────────────────────────────────
# 7. Write outputs
# ─────────────────────────────────────────────────────────────────────
def write_fasta(path, entries):
    with open(path, 'w') as f:
        for h, s in entries:
            f.write(f'>{h}\n')
            for i in range(0, len(s), 60):
                f.write(s[i:i+60] + '\n')

def write_codons(path, entries):
    with open(path, 'w') as f:
        for h, s in entries:
            f.write(f'>{h}\n')
            codons = [s[i:i+3] for i in range(0, len(s), 3)]
            f.write(' '.join(codons) + '\n')

def write_metadata(path, entries):
    with open(path, 'w') as f:
        f.write('header\tnt_len\taa_len\tgc_pct\n')
        for h, s in entries:
            gc = 100.0 * (s.count('G') + s.count('C')) / max(1, len(s))
            f.write(f'{h}\t{len(s)}\t{len(s)//3}\t{gc:.2f}\n')

def write_proteins(path, entries):
    with open(path, 'w') as f:
        for h, s in entries:
            prot = translate(s)
            # Trim trailing stop
            if prot.endswith('*'):
                prot = prot[:-1]
            f.write(f'>{h}\n')
            for i in range(0, len(prot), 60):
                f.write(prot[i:i+60] + '\n')

# ─────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────
def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else 'acquire_homologs/outputs_expanded'
    os.makedirs(out_dir, exist_ok=True)
    api_key = os.environ.get('NCBI_API_KEY', '')

    # Load parent sequence
    parent_fasta = '/resnick/groups/shapirolab/trao2/ndh2_design/parent_info/parent_enzyme.fasta'
    with open(parent_fasta) as f:
        lines = f.read().splitlines()
    parent_seq = ''.join(l for l in lines if not l.startswith('>'))

    # ── Step 1: BLAST ──
    print(f"[{ts()}] === Step 1: BLAST parent against nr (bacteria) ===")
    rid = submit_blast(parent_seq, max_hits=5000)
    blast_text = poll_blast(rid, max_wait=3600)
    # Save raw BLAST output
    with open(os.path.join(out_dir, 'blast_raw.txt'), 'w') as f:
        f.write(blast_text)
    blast_accs = parse_blast_accessions(blast_text)
    print(f"[{ts()}] BLAST accessions parsed: {len(blast_accs)}")

    # Convert accessions to protein IDs
    print(f"[{ts()}] Converting BLAST accessions to protein IDs...")
    blast_prot_ids = accessions_to_ids(blast_accs, api_key=api_key)
    print(f"[{ts()}] BLAST protein IDs: {len(blast_prot_ids)}")

    # ── Step 2: Broad ESearch ──
    print(f"[{ts()}] === Step 2: Broad ESearch queries ===")
    esearch_ids = esearch_protein_broad(api_key=api_key)
    print(f"[{ts()}] ESearch protein IDs: {len(esearch_ids)}")

    # ── Combine ──
    all_prot_ids = list(set(blast_prot_ids) | esearch_ids)
    print(f"[{ts()}] Combined unique protein IDs: {len(all_prot_ids)}")

    # ── Step 3: Link to nuccore ──
    print(f"[{ts()}] === Step 3: Link protein -> nuccore ===")
    nuccore_ids = elink_to_nuccore(all_prot_ids, api_key=api_key)
    print(f"[{ts()}] Nuccore IDs: {len(nuccore_ids)}")

    if not nuccore_ids:
        print('ERROR: No nuccore IDs found')
        sys.exit(1)

    # Write nuccore IDs for transparency/resume debugging
    with open(os.path.join(out_dir, 'nuccore_ids.txt'), 'w') as f:
        f.write('\n'.join(nuccore_ids) + '\n')

    # ── Step 4: Fetch CDS DNA ──
    print(f"[{ts()}] === Step 4: Fetch CDS DNA ===")
    fasta_text = efetch_cds_na(nuccore_ids, out_dir=out_dir, api_key=api_key)
    entries = parse_fasta(fasta_text)
    print(f"[{ts()}] Total CDS entries fetched: {len(entries)}")

    # ── Step 5: Filter ──
    print(f"[{ts()}] === Step 5: Filter ===")
    # First pass: NDH2-like header OR keep all from BLAST (more permissive)
    ndh2_entries = [(h, s) for h, s in entries if is_ndh2_like(h) or not is_excluded(h)]
    print(f"[{ts()}] After header pre-filter: {len(ndh2_entries)}")

    filtered = filter_sequences(ndh2_entries, min_aa=330, max_aa=500)
    print(f"[{ts()}] After length/quality filter (330-500 aa): {len(filtered)}")

    deduped = dedup_sequences(filtered)
    print(f"[{ts()}] After exact dedup: {len(deduped)}")

    # ── Step 6: Also merge previous results ──
    prev_fasta = 'acquire_homologs/outputs/homologs_filtered.fasta'
    if os.path.exists(prev_fasta):
        prev = parse_fasta(open(prev_fasta).read())
        prev_filtered = filter_sequences(prev, min_aa=330, max_aa=500)
        combined = deduped + prev_filtered
        deduped = dedup_sequences(combined)
        print(f"[{ts()}] After merging previous results: {len(deduped)}")

    # ── Step 7: Write outputs ──
    print(f"[{ts()}] === Step 6: Writing outputs ===")
    write_fasta(os.path.join(out_dir, 'ndh2_cds_expanded.fasta'), deduped)
    write_codons(os.path.join(out_dir, 'ndh2_cds_expanded_codons.txt'), deduped)
    write_metadata(os.path.join(out_dir, 'ndh2_cds_expanded_metadata.tsv'), deduped)
    write_proteins(os.path.join(out_dir, 'ndh2_proteins_expanded.fasta'), deduped)

    # Summary
    lengths = [len(s)//3 for _, s in deduped]
    summary = [
        f"Total unique CDS sequences: {len(deduped)}",
        f"AA length: mean {sum(lengths)/len(lengths):.1f}, min {min(lengths)}, max {max(lengths)}",
        f"Outputs in: {out_dir}/",
    ]
    print(f"\n[{ts()}] === DONE ===")
    for l in summary:
        print(f"  {l}")
    with open(os.path.join(out_dir, 'expansion_summary.txt'), 'w') as f:
        f.write('\n'.join(summary) + '\n')

if __name__ == '__main__':
    main()
