# DeepPocket Adapter

DeepPocket is included as a learned rescoring and segmentation layer on top of
fpocket candidate pockets. The adapter runs the upstream `predict.py`, keeps raw
DeepPocket/fpocket outputs, and projects DeepPocket pocket confidence onto the
benchmark residue score TSV format.

Sources:

- Project: https://github.com/devalab/DeepPocket
- Paper: https://doi.org/10.1021/acs.jcim.1c00799

## Setup

DeepPocket requires fpocket, PyTorch, libmolgrid/molgrid, Biopython, scikit-image
and the trained classifier and segmentation checkpoints. The upstream README
points to downloadable model files:

- `first_model_fold1_best_test_auc_85001.pth.tar`
- `seg0_best_test_IOU_91.pth.tar`

Example layout:

```bash
mkdir -p ~/.cache/deeppocket/models
# place the upstream checkpoint files here:
# ~/.cache/deeppocket/models/first_model_fold1_best_test_auc_85001.pth.tar
# ~/.cache/deeppocket/models/seg0_best_test_IOU_91.pth.tar
```

On this aarch64 host, the local DeepPocket run uses `envs/deeppocket.Dockerfile`
to build a `linux/amd64` CPU image. The image follows the upstream Docker
requirements, builds fpocket from source, installs the PyPI `molgrid` wheel, and
patches DeepPocket's CUDA-only calls to run on CPU.

```bash
sudo -n docker build --platform linux/amd64 \
  -f envs/deeppocket.Dockerfile \
  -t binding-deeppocket:latest .
```

## Smoke Run

Run the 20-target smoke manifest:

```bash
python scripts/run_deeppocket.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --out-dir results/deeppocket/p2rank_coach420_mlig_smoke20 \
  --docker-image binding-deeppocket:latest \
  --docker-sudo \
  --class-checkpoint ~/.cache/deeppocket/models/first_model_fold1_best_test_auc_85001.pth.tar \
  --seg-checkpoint ~/.cache/deeppocket/models/seg0_best_test_IOU_91.pth.tar \
  --rank 0 \
  --timeout 900 \
  --resume
```

For a one-target sanity check:

```bash
python scripts/run_deeppocket.py \
  --target-id 1a26A \
  --docker-image binding-deeppocket:latest \
  --docker-sudo \
  --class-checkpoint ~/.cache/deeppocket/models/first_model_fold1_best_test_auc_85001.pth.tar \
  --seg-checkpoint ~/.cache/deeppocket/models/seg0_best_test_IOU_91.pth.tar \
  --rank 0 \
  --timeout 900
```

## Score Normalization

DeepPocket ranks fpocket candidates and writes:

- `raw/<target_id>/<target_id>_nowat_out/pockets/bary_centers_ranked.types`
- `raw/<target_id>/<target_id>_nowat_out/pockets/bary_centers_confidence.txt`
- `raw/<target_id>/<target_id>_nowat_out/pockets/pocket*_atm.pdb`
- optional segmented `raw/<target_id>/<target_id>_nowat_pocket*.pdb`

By default, the adapter uses `--projection-source ranked-candidates`: each
residue in a ranked fpocket candidate receives that candidate's DeepPocket
classifier confidence, and residues in multiple pockets keep their maximum
confidence.

Use `--projection-source segmented-pockets` to project only segmented pocket PDB
files. In that mode, pockets are ranked by file order because the upstream
segmented PDB files do not carry per-residue probabilities.

The smoke20 run used `--rank 0` and the default `ranked-candidates` projection,
so the segmentation model file is provided for upstream argument validation but
segmentation is not executed. If the upstream run completes but no residues are
assigned to any candidate, the adapter writes an all-zero score file.

The normalized score TSV contains `residue_id`, `chain_id`, `residue_number`,
`probability`, `is_binding`, `deeppocket_score`, and `pocket`.

## Residue-Level Evaluation

```bash
python scripts/evaluate_residue_scores.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --scores-dir results/deeppocket/p2rank_coach420_mlig_smoke20/scores \
  --out-dir results/deeppocket/p2rank_coach420_mlig_smoke20/evaluation_residue \
  --distance-cutoff 4.0 \
  --threshold 0.5
```

Smoke20 Docker rank-only result:

- Targets: 20
- Residues: 5025, binding residues: 322
- Pooled AUROC: 0.8888
- Pooled average precision: 0.4139
- Mean AUROC: 0.8783
- Mean average precision: 0.5534
- Mean F1 at threshold 0.5: 0.4869
