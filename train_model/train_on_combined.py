#!/usr/bin/env python3
import os
import sys
import json
import math
import random
import torch
from typing import List
from torch.utils.data import Dataset
from transformers import GPT2LMHeadModel, Trainer, TrainingArguments, GPT2Config

# Ensure repo root is on path
sys.path.insert(0, '/resnick/groups/shapirolab/trao2/ndh2_design')


def load_vocab(vocab_path: str):
    with open(vocab_path, 'r') as f:
        return json.load(f)


def seq_to_ids(seq: str, vocab: dict, cls_id: int, sep_id: int, unk_id: int, max_len: int):
    tokens = seq.strip().split()
    token_ids: List[int] = []
    for tok in tokens:
        tid = vocab.get(tok, unk_id)
        # Guard invalid ids
        if isinstance(tid, int):
            token_ids.append(tid)
        else:
            token_ids.append(unk_id)
    input_ids = [cls_id] + token_ids + [sep_id]
    if len(input_ids) > max_len:
        input_ids = input_ids[:max_len-1] + [sep_id]
    return input_ids


class CodonDataset(Dataset):
    def __init__(self, path: str, vocab: dict, max_len: int = 512):
        self.vocab = vocab
        self.max_len = max_len
        self.unk_id = vocab.get('[UNK]', 0)
        self.cls_id = vocab.get('[CLS]', 1)
        self.sep_id = vocab.get('[SEP]', 2)
        self.pad_id = vocab.get('[PAD]', 3)
        with open(path, 'r') as f:
            self.lines = [ln.strip() for ln in f if ln.strip()]

    def __len__(self):
        return len(self.lines)

    def __getitem__(self, idx):
        input_ids = seq_to_ids(self.lines[idx], self.vocab, self.cls_id, self.sep_id, self.unk_id, self.max_len)
        return {
            'input_ids': torch.tensor(input_ids, dtype=torch.long),
            'attention_mask': torch.tensor([1] * len(input_ids), dtype=torch.long),
            'labels': torch.tensor(input_ids, dtype=torch.long),
        }


def collate_pad(batch, pad_id: int):
    max_len = max(len(item['input_ids']) for item in batch)
    input_ids, attention_mask, labels = [], [], []
    for item in batch:
        cur = item['input_ids']
        pad_n = max_len - len(cur)
        input_ids.append(torch.nn.functional.pad(cur, (0, pad_n), value=pad_id))
        attention_mask.append(torch.tensor([1] * len(cur) + [0] * pad_n, dtype=torch.long))
        labels.append(torch.nn.functional.pad(item['labels'], (0, pad_n), value=pad_id))
    return {
        'input_ids': torch.stack(input_ids),
        'attention_mask': torch.stack(attention_mask),
        'labels': torch.stack(labels),
    }


def main():
    # Args: [train_txt] [val_txt] [model_dir] [out_dir]
    train_txt = sys.argv[1] if len(sys.argv) > 1 else 'acquire_homologs/combined/combined_train.txt'
    val_txt = sys.argv[2] if len(sys.argv) > 2 else 'acquire_homologs/combined/combined_val.txt'
    model_dir = sys.argv[3] if len(sys.argv) > 3 else 'genslm_25M_local'
    out_dir = sys.argv[4] if len(sys.argv) > 4 else 'train_model/checkpoints_combined'

    os.makedirs(out_dir, exist_ok=True)

    random.seed(42)
    torch.manual_seed(42)

    # Load model and vocab
    print(f'[train] Loading model from {model_dir} ...')
    try:
        model = GPT2LMHeadModel.from_pretrained(model_dir, local_files_only=True)
        print('[train] Loaded model via from_pretrained() using default weights file')
    except Exception as e:
        print(f"[train] Primary load failed: {e}")
        # Fallback: load config then explicit state dict if available
        cfg_path = os.path.join(model_dir, 'config.json')
        fixed_path = os.path.join(model_dir, 'pytorch_model_fixed.bin')
        if os.path.isfile(cfg_path) and os.path.isfile(fixed_path):
            print('[train] Falling back to config + pytorch_model_fixed.bin ...')
            config = GPT2Config.from_pretrained(model_dir, local_files_only=True)
            model = GPT2LMHeadModel(config)
            state = torch.load(fixed_path, map_location='cpu')
            missing, unexpected = model.load_state_dict(state, strict=False)
            print(f'[train] Loaded fixed weights. Missing keys: {len(missing)}, Unexpected keys: {len(unexpected)}')
        else:
            raise
    vocab = load_vocab(os.path.join(model_dir, 'vocab.json'))
    pad_id = vocab.get('[PAD]', 3)

    # Datasets
    max_len = 512
    train_ds = CodonDataset(train_txt, vocab, max_len=max_len)
    val_ds = CodonDataset(val_txt, vocab, max_len=max_len)

    print(f'[train] Train N={len(train_ds)} | Val N={len(val_ds)} | max_len={max_len}')

    # Collator
    def _collate(batch):
        return collate_pad(batch, pad_id)

    # Training args
    args = TrainingArguments(
        output_dir=out_dir,
        num_train_epochs=10,
        per_device_train_batch_size=2,
        gradient_accumulation_steps=4,
        learning_rate=5e-5,
        fp16=torch.cuda.is_available(),
        logging_steps=50,
        save_steps=500,
        evaluation_strategy='steps',
        eval_steps=500,
        save_total_limit=3,
        report_to='none',
        overwrite_output_dir=True,
    )

    trainer = Trainer(
        model=model,
        args=args,
        train_dataset=train_ds,
        eval_dataset=val_ds,
        data_collator=_collate,
    )

    print('[train] Starting fine-tuning on combined dataset ...')
    trainer.train()

    final_dir = os.path.join('train_model', 'model_combined_final')
    os.makedirs(final_dir, exist_ok=True)
    trainer.save_model(final_dir)
    print(f'[train] Training complete. Final model saved to {final_dir}')


if __name__ == '__main__':
    main()
