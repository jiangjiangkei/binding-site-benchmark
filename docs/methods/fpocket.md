# fpocket Adapter

fpocket is included as a fast geometry-based baseline for pocket detection. The
adapter runs the fpocket CLI, keeps raw pocket outputs, and projects pocket
scores onto residues so the existing residue-level evaluation and ROC/PR plots
can be reused.

Sources:

- Current project: https://github.com/Discngine/fpocket
- fpocket 1/2 archive: https://sourceforge.net/projects/fpocket/files/

Local smoke run used:

- `fpocket2.tar.gz` from SourceForge
- Executable: `fpocket`
- Input format: PDB

## Setup

On systems with a working fpocket package or binary, use that executable
directly with `--fpocket-bin`.

On this aarch64 machine, the current fpocket 4.x source release needs a molfile
plugin library for the local architecture. For the PDB-only smoke benchmark, the
older fpocket2 archive is enough and avoids that dependency:

```bash
mkdir -p ~/.cache/fpocket
cd ~/.cache/fpocket
curl -L -o fpocket2.tar.gz \
  https://downloads.sourceforge.net/project/fpocket/fpocket2.tar.gz
mkdir -p fpocket2
tar -xzf fpocket2.tar.gz -C fpocket2 --strip-components=1
cd fpocket2
make
```

With newer GCC versions, fpocket2 may need local build-only compatibility
patches for implicit function declarations, link order, and output-path string
formatting. Those patches are not part of this benchmark repository; the
benchmark adapter only requires a callable `fpocket` executable.

## Smoke Run

Run the 20-target smoke manifest:

```bash
python scripts/run_fpocket.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --out-dir results/fpocket/p2rank_coach420_mlig_smoke20 \
  --fpocket-bin ~/.cache/fpocket/fpocket2/bin/fpocket \
  --resume
```

For a one-target sanity check:

```bash
python scripts/run_fpocket.py \
  --target-id 1a26A \
  --fpocket-bin ~/.cache/fpocket/fpocket2/bin/fpocket
```

## Score Normalization

fpocket writes pocket-level scores and per-pocket contacted atoms. The adapter:

1. reads `raw/<target_id>/<target_id>_out/<target_id>_info.txt`;
2. reads `raw/<target_id>/<target_id>_out/pockets/pocket*_atm.pdb`;
3. assigns each residue the highest-scoring pocket it appears in;
4. writes `scores/<target_id>.scores.tsv`.

By default, `probability` is `fpocket_score / max_fpocket_score` within each
target, so residue ranking follows fpocket's native pocket score. Pass
`--probability-source druggability` to project fpocket's druggability score
instead.

The normalized score TSV contains `residue_id`, `chain_id`, `residue_number`,
`probability`, `is_binding`, `fpocket_score`, `druggability_score`, and
`pocket`.

## Residue-Level Evaluation

Evaluate against ligand-contact residues:

```bash
python scripts/evaluate_residue_scores.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --scores-dir results/fpocket/p2rank_coach420_mlig_smoke20/scores \
  --out-dir results/fpocket/p2rank_coach420_mlig_smoke20/evaluation_residue \
  --distance-cutoff 4.0 \
  --threshold 0.5
```

Plot fpocket together with ProtCross and P2Rank:

```bash
python scripts/plot_residue_curves.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --out-dir results/plots/p2rank_coach420_mlig_smoke20/evaluation_residue \
  --method ProtCross=results/protcross/p2rank_coach420_mlig_smoke20/scores \
  --method P2Rank=results/p2rank/p2rank_coach420_mlig_smoke20/scores \
  --method fpocket=results/fpocket/p2rank_coach420_mlig_smoke20/scores
```

This writes the comparison `roc_curve.*`, `pr_curve.*`, and
`curve_metrics.csv` files.
