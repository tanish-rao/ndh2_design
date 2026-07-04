# Fine-tune GenSLM on Combined NDH2 Dataset

This folder contains scripts to fine-tune the local GenSLM model on the combined NDH2 codon-token dataset (Strict+HMM ∪ HMM-filtered balanced set).

## Inputs
- acquire_homologs/combined/combined_train.txt
- acquire_homologs/combined/combined_val.txt
- genslm_25M_local/ (pretrained weights + vocab.json)

## Train (Slurm)
```
sbatch train_model/run_train_combined.sbatch
```
Logs are written to `train_model/logs/`.

## Train (interactive)
```
python -u train_model/train_on_combined.py \
  acquire_homologs/combined/combined_train.txt \
  acquire_homologs/combined/combined_val.txt \
  genslm_25M_local \
  train_model/checkpoints_combined
```

## Outputs
- train_model/checkpoints_combined/ — checkpoints during training
- train_model/model_combined_final/ — final fine-tuned model

## Notes
- The trainer pads with the `[PAD]` id from the provided vocab and trains with Causal LM objective.
- Default hyperparameters: 10 epochs, batch size 2, grad-accum 4, LR 5e-5; adjust in `train_on_combined.py`.
