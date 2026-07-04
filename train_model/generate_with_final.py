#!/usr/bin/env python3
import os
import sys
import json
import math
import time
import argparse
import torch
from transformers import GPT2LMHeadModel


def load_vocab(path):
    with open(path, 'r') as f:
        return json.load(f)


def invert_vocab(vocab):
    return {int(v): k for k, v in vocab.items() if isinstance(v, int) or str(v).isdigit()}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--model_dir', default='train_model/model_combined_final')
    p.add_argument('--vocab_path', default='genslm_25M_local/vocab.json')
    p.add_argument('--out_dir', default='train_model/generation')
    p.add_argument('--num_sequences', type=int, default=100)
    p.add_argument('--batch_size', type=int, default=8)
    p.add_argument('--max_codons', type=int, default=512)
    p.add_argument('--temperature', type=float, default=0.9)
    p.add_argument('--top_p', type=float, default=0.95)
    p.add_argument('--top_k', type=int, default=0)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--prompt_codons', type=str, default='')
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    torch.manual_seed(args.seed)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = GPT2LMHeadModel.from_pretrained(args.model_dir, local_files_only=True).to(device)
    vocab = load_vocab(args.vocab_path)
    id2tok = invert_vocab(vocab)

    cls_id = vocab.get('[CLS]', 1)
    sep_id = vocab.get('[SEP]', 2)
    pad_id = vocab.get('[PAD]', 3)

    if args.prompt_codons.strip():
        toks = args.prompt_codons.strip().split()
        input_ids = [vocab.get(t, vocab.get('[UNK]', 0)) for t in toks]
        input_ids = [cls_id] + input_ids
    else:
        input_ids = [cls_id]

    input_tensor = torch.tensor([input_ids], dtype=torch.long, device=device)

    ts = time.strftime('%Y%m%d_%H%M%S')
    out_path = os.path.join(args.out_dir, f'generated_codons_{ts}.txt')

    model.eval()
    results = []
    remaining = int(args.num_sequences)
    batch = max(1, int(args.batch_size))
    while remaining > 0:
        cur_b = min(batch, remaining)
        with torch.no_grad():
            gen = model.generate(
                input_ids=input_tensor,
                do_sample=True,
                temperature=args.temperature,
                top_p=args.top_p,
                top_k=args.top_k if args.top_k > 0 else None,
                max_new_tokens=args.max_codons,
                eos_token_id=sep_id,
                pad_token_id=pad_id,
                bos_token_id=cls_id,
                num_return_sequences=cur_b,
            )
        # Parse each returned sequence
        for i in range(cur_b):
            seq_ids = gen[i].tolist()
            if seq_ids and seq_ids[0] == cls_id:
                seq_ids = seq_ids[1:]
            if sep_id in seq_ids:
                seq_ids = seq_ids[:seq_ids.index(sep_id)]
            toks = [id2tok.get(tid, '[UNK]') for tid in seq_ids]
            toks = [t for t in toks if t not in ('[UNK]', '[CLS]', '[SEP]', '[PAD]', '[MASK]')]
            results.append(' '.join(toks))
        remaining -= cur_b

    with open(out_path, 'w') as f:
        for r in results:
            if r.strip():
                f.write(r.strip() + '\n')

    print(f'Wrote {len(results)} sequences to {out_path}')


if __name__ == '__main__':
    main()
