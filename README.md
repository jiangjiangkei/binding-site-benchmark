# Binding Site Benchmark

Workspace for building a reproducible protein-small molecule binding-site
prediction benchmark.

## Contents

- `docs/internal-binding-site-benchmark-blueprint.md`: internal v1 benchmark
  blueprint.
- `docs/methods/protcross.md`: ProtCross adapter notes, environment setup,
  smoke-run commands, and residue-level evaluation commands.
- `envs/protcross.yml`: dedicated Conda environment for ProtCross runs.
- `requirements/protcross.txt`: pip requirements for the ProtCross adapter.
- `scripts/run_protcross.py`: batch runner around the ProtCross CLI.
- `scripts/evaluate_residue_scores.py`: residue-level metric evaluation.
- `scripts/plot_residue_curves.py`: pooled ROC and precision-recall plots.

## Data Policy

Benchmark data snapshots and generated outputs stay local. The repository
explicitly ignores raw structures, split files, manifests, runtime assets, and
`results/` outputs in `.gitignore`.

Expected local paths include:

- `data/raw/`
- `data/manifests/`
- `data/splits/`
- `data/assets/`
- `results/`
