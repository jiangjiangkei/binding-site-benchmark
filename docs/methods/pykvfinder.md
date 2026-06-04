# pyKVFinder Adapter

pyKVFinder is included as an open-source geometric cavity-detection baseline.
The adapter runs pyKVFinder on a protein-only PDB, keeps raw cavity summaries,
and projects cavity residue assignments onto benchmark residue score TSVs.

Sources:

- Project: https://github.com/LBC-LNBio/pyKVFinder
- PyPI package: https://pypi.org/project/pyKVFinder/

## Setup

pyKVFinder is available from PyPI and installs cleanly in a lightweight Python
environment:

```bash
python3 -m venv ~/.cache/pykvfinder/venv
~/.cache/pykvfinder/venv/bin/python -m pip install -U pip
~/.cache/pykvfinder/venv/bin/python -m pip install pyKVFinder==0.9.2
```

## Smoke Run

Run the 20-target smoke manifest:

```bash
python scripts/run_pykvfinder.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --out-dir results/pykvfinder/p2rank_coach420_mlig_smoke20 \
  --python-bin ~/.cache/pykvfinder/venv/bin/python \
  --nthreads 1 \
  --timeout 120 \
  --resume
```

For a one-target sanity check:

```bash
python scripts/run_pykvfinder.py \
  --target-id 1a26A \
  --python-bin ~/.cache/pykvfinder/venv/bin/python \
  --nthreads 1 \
  --timeout 120
```

## Score Normalization

The adapter writes a protein-only PDB before calling pyKVFinder. This keeps
benchmark ligands and waters out of the cavity-detection input.

Raw per-target outputs are:

- `raw/<target_id>/<target_id>_protein.pdb`
- `raw/<target_id>/pykvfinder_results.json`
- `raw/<target_id>/pykvfinder_results.toml`

By default, `--score-mode volume` assigns each cavity `volume / max_volume`
within the same target. Residues assigned to multiple cavities keep their
maximum probability. Additional modes are `depth`, `rank-decay`, and `binary`.
If pyKVFinder completes but detects no cavities, the adapter writes an all-zero
score file for that target.

The normalized score TSV contains `residue_id`, `chain_id`, `residue_number`,
`probability`, `is_binding`, `pykvfinder_score`, and `pocket`.

## Residue-Level Evaluation

```bash
python scripts/evaluate_residue_scores.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --scores-dir results/pykvfinder/p2rank_coach420_mlig_smoke20/scores \
  --out-dir results/pykvfinder/p2rank_coach420_mlig_smoke20/evaluation_residue \
  --distance-cutoff 4.0 \
  --threshold 0.5
```

Smoke20 result:

- Targets: 20
- Residues: 5025, binding residues: 322
- Pooled AUROC: 0.8596
- Pooled average precision: 0.3120
- Mean AUROC: 0.8383
- Mean average precision: 0.3967
- Mean F1 at threshold 0.5: 0.4319
