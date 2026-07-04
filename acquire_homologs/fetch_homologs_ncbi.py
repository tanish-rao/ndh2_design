#!/usr/bin/env python3
import os
import sys
import time
import math
import json
import re
from typing import List, Dict, Tuple
from datetime import datetime

try:
    import requests
    from requests.exceptions import ChunkedEncodingError, ConnectionError, ReadTimeout
except Exception:
    print('ERROR: This script requires the requests package. Install via: pip install requests')
    sys.exit(1)

NCBI_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
UA = {'User-Agent': 'ndh2-homolog-fetcher/1.0'}

def _fetch_json(url: str, params: dict, timeout: int = 60) -> dict:
    delay = 1.0
    for attempt in range(6):
        try:
            r = requests.get(url, params=params, headers=UA, timeout=timeout)
            r.raise_for_status()
            return r.json()
        except (ChunkedEncodingError, ConnectionError, ReadTimeout):
            if attempt == 5:
                raise
            time.sleep(delay)
            delay *= 2

def _fetch_text(url: str, params: dict, timeout: int = 120) -> str:
    delay = 1.0
    for attempt in range(6):
        try:
            r = requests.get(url, params=params, headers=UA, timeout=timeout)
            r.raise_for_status()
            return r.text
        except (ChunkedEncodingError, ConnectionError, ReadTimeout):
            if attempt == 5:
                raise
            time.sleep(delay)
            delay *= 2

def esearch_protein(term: str, retmax: int = 500, api_key: str = '') -> List[str]:
    params = {
        'db': 'protein',
        'retmode': 'json',
        'retmax': str(retmax),
        'term': term,
    }
    if api_key:
        params['api_key'] = api_key
    data = _fetch_json(f'{NCBI_BASE}/esearch.fcgi', params, timeout=60)
    return data.get('esearchresult', {}).get('idlist', [])

def elink_protein_to_nuccore(prot_ids: List[str], api_key: str = '') -> List[str]:
    nuccore_ids = set()
    B = 200
    for i in range(0, len(prot_ids), B):
        batch = prot_ids[i:i+B]
        params = {
            'dbfrom': 'protein',
            'db': 'nuccore',
            'id': ','.join(batch),
            'retmode': 'json',
        }
        if api_key:
            params['api_key'] = api_key
        data = _fetch_json(f'{NCBI_BASE}/elink.fcgi', params, timeout=60)
        for linkset in data.get('linksets', []):
            for linksetdb in linkset.get('linksetdbs', []) or []:
                for link in linksetdb.get('links', []) or []:
                    nuccore_ids.add(str(link))
        time.sleep(0.34)
    return list(nuccore_ids)

def efetch_fasta_cds_na(nuccore_ids: List[str], api_key: str = '') -> str:
    fasta = ''
    B = 25
    for i in range(0, len(nuccore_ids), B):
        batch = nuccore_ids[i:i+B]
        params = {
            'db': 'nuccore',
            'rettype': 'fasta_cds_na',
            'retmode': 'text',
            'id': ','.join(batch),
        }
        if api_key:
            params['api_key'] = api_key
        fasta += _fetch_text(f'{NCBI_BASE}/efetch.fcgi', params, timeout=180)
        time.sleep(0.34)
    return fasta

def parse_fasta(text: str) -> List[Tuple[str, str]]:
    entries = []
    header = None
    seq_lines = []
    for line in text.splitlines():
        if not line:
            continue
        if line.startswith('>'):
            if header is not None:
                entries.append((header, ''.join(seq_lines).replace(' ', '').replace('\t', '')))
            header = line[1:].strip()
            seq_lines = []
        else:
            seq_lines.append(line.strip())
    if header is not None:
        entries.append((header, ''.join(seq_lines).replace(' ', '').replace('\t', '')))
    return entries

def header_is_ndh2(h: str) -> bool:
    hs = h.lower()
    if 'hypothetical' in hs:
        return False
    if 'partial cds' in hs:
        return False
    if 'ndh2' in hs or 'ndh-2' in hs:
        return True
    if 'nadh dehydrogenase (quinone)' in hs:
        return True
    if 'type ii nadh dehydrogenase' in hs:
        return True
    return False

def filter_length(seqs: List[Tuple[str,str]], min_nt=1000, max_nt=1500) -> List[Tuple[str,str]]:
    out = []
    for h, s in seqs:
        if len(s) % 3 != 0:
            continue
        if min_nt <= len(s) <= max_nt:
            out.append((h, s))
    return out

def dedup_exact(seqs: List[Tuple[str,str]]) -> List[Tuple[str,str]]:
    seen = {}
    out = []
    for h, s in seqs:
        if s in seen:
            continue
        seen[s] = True
        out.append((h, s))
    return out

def write_fasta(path: str, seqs: List[Tuple[str,str]]):
    with open(path, 'w') as f:
        for h, s in seqs:
            f.write(f'>{h}\n')
            for i in range(0, len(s), 60):
                f.write(s[i:i+60] + '\n')

def write_metadata(path: str, seqs: List[Tuple[str,str]]):
    with open(path, 'w') as f:
        f.write('header\tlength_nt\tlength_aa\n')
        for h, s in seqs:
            f.write(f'{h}\t{len(s)}\t{len(s)//3}\n')

def write_codon_triplets(path: str, seqs: List[Tuple[str,str]]):
    with open(path, 'w') as f:
        for h, s in seqs:
            f.write(f'>{h}\n')
            codons = [s[i:i+3] for i in range(0, len(s), 3)]
            f.write(' '.join(codons) + '\n')

def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else 'acquire_homologs/outputs'
    os.makedirs(out_dir, exist_ok=True)
    api_key = os.environ.get('NCBI_API_KEY', '')
    term = '((ndh2[Gene Name]) OR (ndh-2[All Fields]) OR ("NADH dehydrogenase (quinone)"[Title]) OR ("type II NADH dehydrogenase"[Title])) AND bacteria[Organism] NOT hypothetical'
    print(f"[{datetime.now().isoformat(timespec='seconds')}] ESearch proteins...")
    prot_ids = esearch_protein(term, retmax=1000, api_key=api_key)
    print(f"[{datetime.now().isoformat(timespec='seconds')}] Protein IDs: {len(prot_ids)}")
    if not prot_ids:
        print('No protein IDs found for query')
        sys.exit(1)
    print(f"[{datetime.now().isoformat(timespec='seconds')}] Link protein->nuccore...")
    nuccore_ids = elink_protein_to_nuccore(prot_ids, api_key=api_key)
    print(f"[{datetime.now().isoformat(timespec='seconds')}] Nuccore IDs: {len(nuccore_ids)}")
    if not nuccore_ids:
        print('No nuccore IDs linked from protein IDs')
        sys.exit(1)
    print(f"[{datetime.now().isoformat(timespec='seconds')}] EFetch fasta_cds_na...")
    fasta_text = efetch_fasta_cds_na(nuccore_ids, api_key=api_key)
    entries = parse_fasta(fasta_text)
    print(f"[{datetime.now().isoformat(timespec='seconds')}] Entries fetched: {len(entries)}")
    entries = [(h, s) for h, s in entries if header_is_ndh2(h)]
    print(f"[{datetime.now().isoformat(timespec='seconds')}] After NDH2 header filter: {len(entries)}")
    raw_path = os.path.join(out_dir, 'homologs_raw.fasta')
    write_fasta(raw_path, entries)
    filt = filter_length(entries, 1000, 2000)
    print(f"[{datetime.now().isoformat(timespec='seconds')}] After length filter (1000-2000 nt): {len(filt)}")
    dedup = dedup_exact(filt)
    print(f"[{datetime.now().isoformat(timespec='seconds')}] After exact dedup: {len(dedup)}")
    filt_path = os.path.join(out_dir, 'homologs_filtered.fasta')
    write_fasta(filt_path, dedup)
    meta_path = os.path.join(out_dir, 'homologs_metadata.tsv')
    write_metadata(meta_path, dedup)
    codon_path = os.path.join(out_dir, 'homologs_filtered_codons.txt')
    write_codon_triplets(codon_path, dedup)
    print(f'raw: {raw_path}')
    print(f'filtered: {filt_path}')
    print(f'metadata: {meta_path}')
    print(f'codons: {codon_path}')

if __name__ == '__main__':
    main()
