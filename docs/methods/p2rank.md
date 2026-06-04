# P2Rank Adapter

P2Rank is included as a baseline method for residue-level binding-site
prediction by normalizing its per-residue output into the benchmark score
schema.

Source: https://github.com/rdk/p2rank

Pinned local run:

- P2Rank release: `2.5.1`
- Executable: `prank`
- Runtime: Java 17 or later

## Setup

P2Rank is distributed as a portable binary package. One local setup option:

```bash
mkdir -p ~/.cache/p2rank
cd ~/.cache/p2rank
curl -L -o p2rank_2.5.1.tar.gz \
  https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz
tar -xzf p2rank_2.5.1.tar.gz
```

Use the bundled executable directly:

```bash
~/.cache/p2rank/p2rank_2.5.1/prank help
```

## Smoke Run

Run the 20-target smoke manifest:

```bash
python scripts/run_p2rank.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --out-dir results/p2rank/p2rank_coach420_mlig_smoke20 \
  --p2rank-bin ~/.cache/p2rank/p2rank_2.5.1/prank \
  --threads 4 \
  --resume
```

For a one-target sanity check:

```bash
python scripts/run_p2rank.py \
  --target-id 1a26A \
  --p2rank-bin ~/.cache/p2rank/p2rank_2.5.1/prank \
  --threads 1
```

## Outputs

The adapter writes:

- `raw/<input_file>_predictions.csv`: P2Rank pocket-level predictions.
- `raw/<input_file>_residues.csv`: P2Rank residue-level scores.
- `scores/<target_id>.scores.tsv`: normalized benchmark residue scores with
  `residue_id`, `chain_id`, `residue_number`, `probability`, `is_binding`,
  `p2rank_score`, and `pocket`.
- `run_summary.csv`: per-target status, normalized residue counts, command,
  and log paths.

When a manifest row has `chain_id`, normalization keeps only that chain so
P2Rank uses the same residue universe as the other residue-level adapters.

## Residue-Level Evaluation

Evaluate against ligand-contact residues:

```bash
python scripts/evaluate_residue_scores.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --scores-dir results/p2rank/p2rank_coach420_mlig_smoke20/scores \
  --out-dir results/p2rank/p2rank_coach420_mlig_smoke20/evaluation_residue \
  --distance-cutoff 4.0 \
  --threshold 0.5
```

Plot P2Rank together with ProtCross on pooled ROC and PR curves:

```bash
python scripts/plot_residue_curves.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --out-dir results/plots/p2rank_coach420_mlig_smoke20/evaluation_residue \
  --method ProtCross=results/protcross/p2rank_coach420_mlig_smoke20/scores \
  --method P2Rank=results/p2rank/p2rank_coach420_mlig_smoke20/scores
```

This writes `roc_curve.svg`, `roc_curve.png`, `pr_curve.svg`, `pr_curve.png`,
and `curve_metrics.csv` into the comparison output directory.
