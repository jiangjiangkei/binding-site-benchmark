# ProtCross Adapter

ProtCross is included as an experimental benchmark method for residue-level
binding-site prediction.

Source: https://github.com/GeraltZeroZhong/ProtCross

Pinned package:

- `protcross==0.1.2`
- `esm==3.2.1.post1`
- Git tag: `v0.1.2`
- Release checkpoint: `protcross-0.1.2-binding-moad-final.ckpt`
- PCA reducer: `pca_esmc_128_binding_moad_0.1.2.pkl`
- ESM-C weights: `esmc_600m_2024_12_v0.pth`

## Environment

Create the dedicated environment:

```bash
conda env create -f envs/protcross.yml
conda activate binding-protcross
```

If the environment already exists:

```bash
conda env update -f envs/protcross.yml --prune
conda activate binding-protcross
```

## Assets

Install runtime assets into the default ProtCross cache:

```bash
protcross setup-assets
```

The default asset directory is `~/.cache/protcross/assets/v0.1.2`. The ESM-C
weights are about 2.3 GB and are distributed from Hugging Face under
EvolutionaryScale's model terms. On this machine the full asset set is already
installed there. To reuse another local ESM-C weights file:

```bash
protcross setup-assets --skip-esm
export PROTCROSS_ESM_WEIGHTS=/absolute/path/to/esmc_600m_2024_12_v0.pth
```

## Smoke Run

Run the 20-target smoke manifest:

```bash
python scripts/run_protcross.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --out-dir results/protcross/p2rank_coach420_mlig_smoke20 \
  --threshold 0.5 \
  --device auto
```

For a one-target sanity check:

```bash
python scripts/run_protcross.py \
  --target-id 1a26A \
  --limit 1 \
  --device cpu
```

## Outputs

The adapter writes:

- `pdb/<target_id>.protcross.pdb`: input PDB with ProtCross probabilities in
  the B-factor column.
- `scores/<target_id>.scores.tsv`: residue-level scores with
  `residue_id`, `chain_id`, `residue_number`, `probability`, and `is_binding`.
- `run_summary.csv`: per-target command, status, runtime, and log paths.

The benchmark should report the final probability threshold when using
ProtCross residue calls.

## Residue-Level Evaluation

Evaluate against ligand-contact residues:

```bash
python scripts/evaluate_residue_scores.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --scores-dir results/protcross/p2rank_coach420_mlig_smoke20/scores \
  --out-dir results/protcross/p2rank_coach420_mlig_smoke20/evaluation_residue \
  --distance-cutoff 4.0 \
  --threshold 0.5
```

The default ground truth marks a residue positive when any heavy atom from the
manifest chain is within `4.0 A` of any heavy atom from the manifest ligand
codes. This is the benchmark's default residue-level protocol.

Plot pooled ROC and PR curves:

```bash
python scripts/plot_residue_curves.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --scores-dir results/protcross/p2rank_coach420_mlig_smoke20/scores \
  --out-dir results/protcross/p2rank_coach420_mlig_smoke20/evaluation_residue \
  --method-name ProtCross
```
