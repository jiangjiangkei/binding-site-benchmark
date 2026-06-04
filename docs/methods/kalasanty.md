# Kalasanty Adapter

Kalasanty is included as a classic 3D U-Net binding-site segmentation baseline.
The adapter calls the upstream `scripts/predict.py`, keeps predicted pocket
files, and maps pocket atoms back to benchmark residues by coordinate matching.

Sources:

- Project: https://gitlab.com/cheminfIBB/kalasanty
- Paper: https://doi.org/10.1038/s41598-020-61860-z

## Setup

Kalasanty was developed against an older Keras/TensorFlow/Open Babel stack. The
upstream repository includes a pretrained model at `data/model_scpdb2017.hdf`.

Example upstream-style setup:

```bash
git clone https://gitlab.com/cheminfIBB/kalasanty.git ~/.cache/kalasanty/kalasanty
cd ~/.cache/kalasanty/kalasanty
conda env create -f environment.yml -n kalasanty_env
conda activate kalasanty_env
pip install .
```

On this aarch64 host, the base Python environment does not include `pybel`,
`tfbio`, TensorFlow 1.x, or the older Keras stack expected by Kalasanty. The
repository includes `envs/kalasanty.Dockerfile`, which builds a `linux/amd64`
Conda environment with the compatible upstream stack and pretrained model.

```bash
sudo -n docker build --platform linux/amd64 \
  -f envs/kalasanty.Dockerfile \
  -t binding-kalasanty:latest .
```

## Smoke Run

Run the 20-target smoke manifest:

```bash
python scripts/run_kalasanty.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --out-dir results/kalasanty/p2rank_coach420_mlig_smoke20 \
  --docker-image binding-kalasanty:latest \
  --docker-sudo \
  --timeout 900 \
  --resume
```

For a one-target sanity check:

```bash
python scripts/run_kalasanty.py \
  --target-id 1a26A \
  --docker-image binding-kalasanty:latest \
  --docker-sudo \
  --timeout 900
```

The adapter writes a protein-only PDB per target before calling Kalasanty. This
keeps benchmark ligands and waters out of the prediction input.

## Score Normalization

Kalasanty writes:

```text
raw/<target_id>/predicted_pockets/<target_id>/pocket*.mol2
raw/<target_id>/predicted_pockets/<target_id>/pockets.cmap
```

The pocket MOL2 files contain predicted pocket atoms but no per-pocket
confidence score. The adapter maps MOL2 atom coordinates back to the original
PDB atoms and assigns matched residues probability `1.0` by default. Pass
`--pocket-score-mode rank-decay` to assign `1 / rank` by pocket file order.
If the upstream run completes but emits only `pockets.cmap` and no pocket MOL2
files for a target, the adapter writes an all-zero score file for that target.

The normalized score TSV contains `residue_id`, `chain_id`, `residue_number`,
`probability`, `is_binding`, `kalasanty_score`, and `pocket`.

## Residue-Level Evaluation

```bash
python scripts/evaluate_residue_scores.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --scores-dir results/kalasanty/p2rank_coach420_mlig_smoke20/scores \
  --out-dir results/kalasanty/p2rank_coach420_mlig_smoke20/evaluation_residue \
  --distance-cutoff 4.0 \
  --threshold 0.5
```

Smoke20 Docker result:

- Targets: 20
- Residues: 5025, binding residues: 322
- Pooled AUROC: 0.8152
- Pooled average precision: 0.2652
- Mean AUROC: 0.7761
- Mean average precision: 0.3115
- Mean F1 at threshold 0.5: 0.3508
