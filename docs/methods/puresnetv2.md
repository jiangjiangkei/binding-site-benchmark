# PUResNetV2 Adapter

PUResNetV2 is included as a sparse-convolution deep learning baseline for
protein-ligand binding-site prediction. The adapter calls the upstream
`puresnet.predict.make_prediction` API and converts predicted pocket PDB files
into residue score TSVs.

Sources:

- Project: https://github.com/jivankandel/PUResNetV2.0
- PyPI package: https://pypi.org/project/puresnet/
- Paper: https://doi.org/10.1186/s13321-024-00865-6

## Setup

The upstream README currently documents `puresnet==0.1`; PyPI also provides
newer `puresnet` releases with the same `make_prediction` API checked for this
adapter. PUResNetV2 requires PyTorch, MinkowskiEngine, Open Babel, and
scikit-learn.

Example upstream-style environment:

```bash
conda create -n puresnetv2 python=3.10 -c conda-forge
conda activate puresnetv2
conda install openblas-devel -c anaconda
conda install pytorch=1.13.0 torchvision=0.14 pytorch-cuda=11.7 -c pytorch -c nvidia
conda install -c "nvidia/label/cuda-11.7.0" cuda-toolkit
export CUDA_HOME="$CONDA_PREFIX"
pip install -U git+https://github.com/NVIDIA/MinkowskiEngine --no-deps
conda install -c conda-forge openbabel
conda install -c anaconda scikit-learn
pip install puresnet==0.1
```

The official Docker image is `jivankandel/puresnet:latest`. It is published as
`linux/amd64` and is large, but it runs on this aarch64 host through Docker
emulation. The adapter's Docker mode also patches `torch.load` to load upstream
CUDA checkpoints on CPU-only hosts.

The upstream package attempts to download residue CIF/SDF templates dynamically.
For reproducible local runs, clone the upstream repository and pass its checked
residue templates explicitly:

```bash
git clone https://github.com/jivankandel/PUResNetV2.0.git /tmp/PUResNetV2.0
```

## Smoke Run

Run the 20-target smoke manifest:

```bash
python scripts/run_puresnetv2.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --out-dir results/puresnetv2/p2rank_coach420_mlig_smoke20 \
  --docker-image jivankandel/puresnet:latest \
  --docker-sudo \
  --device cpu \
  --residues-dir /tmp/PUResNetV2.0/Example/residues \
  --timeout 900 \
  --resume
```

For a one-target sanity check:

```bash
python scripts/run_puresnetv2.py \
  --target-id 1a26A \
  --docker-image jivankandel/puresnet:latest \
  --docker-sudo \
  --device cpu \
  --residues-dir /tmp/PUResNetV2.0/Example/residues \
  --timeout 900
```

## Score Normalization

PUResNetV2 writes clustered predicted pocket residues to PDB files under the
current working directory's `results/<target_id>/` directory. The adapter runs
each target inside `raw/<target_id>/`, so raw outputs land under:

```text
raw/<target_id>/results/<target_id>/*.pdb
```

The upstream PDB output stores atom-level prediction probability in the
occupancy column. The adapter assigns each residue the maximum occupancy among
its predicted atoms and writes non-predicted residues with probability `0`.
If the upstream run completes but emits no clustered pocket PDB files for a
target, the adapter writes an all-zero score file for that target.

The normalized score TSV contains `residue_id`, `chain_id`, `residue_number`,
`probability`, `is_binding`, `puresnetv2_score`, and `pocket`.

## Residue-Level Evaluation

```bash
python scripts/evaluate_residue_scores.py \
  --manifest data/manifests/p2rank_coach420_mlig_smoke20.csv \
  --scores-dir results/puresnetv2/p2rank_coach420_mlig_smoke20/scores \
  --out-dir results/puresnetv2/p2rank_coach420_mlig_smoke20/evaluation_residue \
  --distance-cutoff 4.0 \
  --threshold 0.5
```

Smoke20 CPU Docker result:

- Targets: 20
- Residues: 5025, binding residues: 322
- Pooled AUROC: 0.9194
- Pooled average precision: 0.6931
- Mean AUROC: 0.8956
- Mean average precision: 0.7426
- Mean F1 at threshold 0.5: 0.6867
