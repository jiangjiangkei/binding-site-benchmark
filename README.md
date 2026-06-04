# Binding Site Benchmark

Workspace for building a reproducible protein-small molecule binding-site
prediction benchmark.

## Contents

- `docs/internal-binding-site-benchmark-blueprint.md`: internal v1 benchmark
  blueprint.
- `docs/methods/protcross.md`: ProtCross adapter notes, environment setup,
  smoke-run commands, and residue-level evaluation commands.
- `docs/methods/p2rank.md`: P2Rank adapter notes, local setup, smoke-run
  commands, and comparison-curve commands.
- `docs/methods/fpocket.md`: fpocket adapter notes, local setup, smoke-run
  commands, and comparison-curve commands.
- `docs/methods/deeppocket.md`: DeepPocket adapter notes, setup requirements,
  smoke-run commands, and residue projection details.
- `docs/methods/puresnetv2.md`: PUResNetV2 adapter notes, setup requirements,
  smoke-run commands, and residue projection details.
- `docs/methods/kalasanty.md`: Kalasanty adapter notes, setup requirements,
  smoke-run commands, and residue projection details.
- `docs/methods/pykvfinder.md`: pyKVFinder adapter notes, setup requirements,
  smoke-run commands, and residue projection details.
- `envs/protcross.yml`: dedicated Conda environment for ProtCross runs.
- `envs/deeppocket.Dockerfile`: CPU Docker image for DeepPocket on hosts without
  native `molgrid`/CUDA support.
- `envs/kalasanty.Dockerfile`: Docker image for the older
  Keras/TensorFlow/Open Babel stack required by Kalasanty.
- `requirements/protcross.txt`: pip requirements for the ProtCross adapter.
- `scripts/run_protcross.py`: batch runner around the ProtCross CLI.
- `scripts/run_p2rank.py`: batch runner around the P2Rank CLI that normalizes
  P2Rank residue CSVs into benchmark score TSVs.
- `scripts/run_fpocket.py`: batch runner around the fpocket CLI that projects
  pocket scores onto benchmark residue score TSVs.
- `scripts/run_deeppocket.py`: batch runner around the DeepPocket upstream
  script that projects ranked or segmented pockets onto residue score TSVs.
- `scripts/run_puresnetv2.py`: batch runner around the PUResNetV2 Python API
  that normalizes predicted pocket PDB files into residue score TSVs.
- `scripts/run_kalasanty.py`: batch runner around the Kalasanty upstream script
  that maps predicted pocket MOL2 atoms back to residue score TSVs.
- `scripts/run_pykvfinder.py`: batch runner around pyKVFinder cavity detection
  that normalizes cavity residue assignments into score TSVs.
- `scripts/residue_score_utils.py`: shared helpers for residue projection in
  pocket-style method adapters.
- `scripts/evaluate_residue_scores.py`: residue-level metric evaluation.
- `scripts/plot_residue_curves.py`: pooled ROC and precision-recall plots,
  including repeated `--method NAME=SCORES_DIR` comparison curves.

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

## Comparison Plot

After a method has written `results/<method>/p2rank_coach420_mlig_smoke20/scores`,
it can be added to the shared pooled ROC and PR plot:

```bash
python scripts/plot_residue_curves.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --out-dir results/plots/p2rank_coach420_mlig_smoke20/evaluation_residue \
  --method ProtCross=results/protcross/p2rank_coach420_mlig_smoke20/scores \
  --method P2Rank=results/p2rank/p2rank_coach420_mlig_smoke20/scores \
  --method fpocket=results/fpocket/p2rank_coach420_mlig_smoke20/scores \
  --method DeepPocket=results/deeppocket/p2rank_coach420_mlig_smoke20/scores \
  --method PUResNetV2=results/puresnetv2/p2rank_coach420_mlig_smoke20/scores \
  --method Kalasanty=results/kalasanty/p2rank_coach420_mlig_smoke20/scores \
  --method pyKVFinder=results/pykvfinder/p2rank_coach420_mlig_smoke20/scores
```
