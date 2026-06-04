#!/usr/bin/env python3
"""Evaluate residue-level binding-site scores against ligand-neighbor labels."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO_ROOT / "data/manifests/p2rank_coach420_mlig_smoke20.csv"
DEFAULT_SCORES_DIR = REPO_ROOT / "results/protcross/p2rank_coach420_mlig_smoke20/scores"
DEFAULT_OUT_DIR = REPO_ROOT / "results/protcross/p2rank_coach420_mlig_smoke20/evaluation_residue"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Evaluate residue-level prediction scores using ligand-contact residue labels.",
    )
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST), help="CSV manifest with structure and ligand metadata.")
    parser.add_argument("--scores-dir", default=str(DEFAULT_SCORES_DIR), help="Directory with <target_id>.scores.tsv files.")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR), help="Directory for evaluation tables.")
    parser.add_argument("--score-suffix", default=".scores.tsv", help="Suffix after target_id for score files.")
    parser.add_argument("--distance-cutoff", type=float, default=4.0, help="Ligand-contact cutoff in Angstrom.")
    parser.add_argument("--threshold", type=float, default=0.5, help="Probability threshold for binary calls.")
    parser.add_argument("--target-id", action="append", default=[], help="Only evaluate this target_id. Can be repeated.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    manifest_path = resolve_path(args.manifest, Path.cwd())
    scores_dir = resolve_path(args.scores_dir, Path.cwd())
    out_dir = resolve_path(args.out_dir, Path.cwd())
    out_dir.mkdir(parents=True, exist_ok=True)

    selected = set(args.target_id)
    rows = load_manifest(manifest_path, selected)
    per_target = []
    all_labels: list[int] = []
    all_scores: list[float] = []

    for row in rows:
        target_id = row["target_id"]
        score_path = scores_dir / f"{target_id}{args.score_suffix}"
        structure_path = resolve_path(row["structure_path"], REPO_ROOT)
        ligand_codes = split_codes(row.get("ligand_codes", ""))
        chain_id = row.get("chain_id", "").strip()

        score_records = load_scores(score_path)
        truth = binding_residues_from_pdb(structure_path, ligand_codes, chain_id, args.distance_cutoff)
        labels = [1 if residue_key(record["chain_id"], record["residue_number"]) in truth else 0 for record in score_records]
        scores = [record["probability"] for record in score_records]
        metrics = calculate_metrics(labels, scores, args.threshold)
        metrics.update(
            {
                "target_id": target_id,
                "pdb_id": row.get("pdb_id", ""),
                "chain_id": chain_id,
                "ligand_codes": ",".join(ligand_codes),
                "score_path": str(score_path),
                "structure_path": str(structure_path),
            }
        )
        per_target.append(metrics)
        all_labels.extend(labels)
        all_scores.extend(scores)

    write_csv(out_dir / "per_target_metrics.csv", per_target)
    summary = summarize(per_target, all_labels, all_scores, args.threshold, args.distance_cutoff)
    write_csv(out_dir / "summary_metrics.csv", [summary])
    print_summary(summary, out_dir)
    return 0


def load_manifest(path: Path, selected_target_ids: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"target_id", "chain_id", "ligand_codes", "structure_path"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Manifest missing required columns: {', '.join(sorted(missing))}")
        return [row for row in reader if not selected_target_ids or row["target_id"] in selected_target_ids]


def load_scores(path: Path) -> list[dict[str, str | float]]:
    if not path.exists():
        raise SystemExit(f"Missing score file: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chain_id", "residue_number", "probability"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Score file missing required columns: {path}: {', '.join(sorted(missing))}")
        records = []
        for row in reader:
            records.append(
                {
                    "chain_id": row["chain_id"],
                    "residue_number": row["residue_number"],
                    "probability": float(row["probability"]),
                }
            )
        return records


def binding_residues_from_pdb(
    pdb_path: Path,
    ligand_codes: set[str],
    chain_id: str,
    distance_cutoff: float,
) -> set[tuple[str, str]]:
    protein_atoms = []
    ligand_atoms_same_chain = []
    ligand_atoms_all = []

    with pdb_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            if record not in {"ATOM", "HETATM"}:
                continue
            atom = parse_pdb_atom(line)
            if atom is None or atom["element"] == "H":
                continue

            atom_chain = atom["chain_id"]
            if record == "ATOM" and (not chain_id or atom_chain == chain_id):
                protein_atoms.append(atom)
            elif record == "HETATM" and atom["resname"] in ligand_codes:
                ligand_atoms_all.append(atom)
                if not chain_id or atom_chain == chain_id:
                    ligand_atoms_same_chain.append(atom)

    ligand_atoms = ligand_atoms_same_chain or ligand_atoms_all
    if not ligand_atoms:
        raise SystemExit(f"No ligand atoms found in {pdb_path} for ligand codes: {', '.join(sorted(ligand_codes))}")

    cutoff2 = distance_cutoff * distance_cutoff
    truth = set()
    for protein_atom in protein_atoms:
        residue = residue_key(protein_atom["chain_id"], protein_atom["residue_number"])
        if residue in truth:
            continue
        px, py, pz = protein_atom["xyz"]
        for ligand_atom in ligand_atoms:
            lx, ly, lz = ligand_atom["xyz"]
            if (px - lx) ** 2 + (py - ly) ** 2 + (pz - lz) ** 2 <= cutoff2:
                truth.add(residue)
                break
    return truth


def parse_pdb_atom(line: str) -> dict[str, object] | None:
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None

    atom_name = line[12:16].strip()
    element = line[76:78].strip() or "".join(char for char in atom_name if char.isalpha())[:1]
    residue_number = line[22:26].strip()
    insertion_code = line[26].strip()
    if insertion_code:
        residue_number = f"{residue_number}{insertion_code}"
    return {
        "chain_id": line[21].strip(),
        "resname": line[17:20].strip(),
        "residue_number": residue_number,
        "element": element.upper(),
        "xyz": (x, y, z),
    }


def calculate_metrics(labels: list[int], scores: list[float], threshold: float) -> dict[str, str | int | float]:
    n = len(labels)
    positives = sum(labels)
    predictions = [1 if score >= threshold else 0 for score in scores]
    true_positive = sum(1 for label, pred in zip(labels, predictions) if label and pred)
    false_positive = sum(1 for label, pred in zip(labels, predictions) if not label and pred)
    false_negative = sum(1 for label, pred in zip(labels, predictions) if label and not pred)
    true_negative = sum(1 for label, pred in zip(labels, predictions) if not label and not pred)

    precision = divide(true_positive, true_positive + false_positive)
    recall = divide(true_positive, true_positive + false_negative)
    specificity = divide(true_negative, true_negative + false_positive)
    f1 = divide(2 * precision * recall, precision + recall)
    top_nplus2_recall = recall_at_k(labels, scores, positives + 2)

    return {
        "n_residues": n,
        "n_binding_residues": positives,
        "positive_rate": divide(positives, n),
        "threshold": threshold,
        "predicted_positive": sum(predictions),
        "precision": precision,
        "recall": recall,
        "specificity": specificity,
        "f1": f1,
        "auroc": auroc(labels, scores),
        "average_precision": average_precision(labels, scores),
        "top_nplus2_recall": top_nplus2_recall,
        "max_score": max(scores) if scores else math.nan,
        "mean_score": sum(scores) / n if n else math.nan,
    }


def summarize(
    per_target: list[dict[str, str | int | float]],
    labels: list[int],
    scores: list[float],
    threshold: float,
    distance_cutoff: float,
) -> dict[str, str | int | float]:
    metric_names = [
        "positive_rate",
        "predicted_positive",
        "precision",
        "recall",
        "specificity",
        "f1",
        "auroc",
        "average_precision",
        "top_nplus2_recall",
    ]
    summary: dict[str, str | int | float] = {
        "n_targets": len(per_target),
        "n_residues": len(labels),
        "n_binding_residues": sum(labels),
        "threshold": threshold,
        "distance_cutoff": distance_cutoff,
        "micro_precision": calculate_metrics(labels, scores, threshold)["precision"],
        "micro_recall": calculate_metrics(labels, scores, threshold)["recall"],
        "micro_f1": calculate_metrics(labels, scores, threshold)["f1"],
        "pooled_auroc": auroc(labels, scores),
        "pooled_average_precision": average_precision(labels, scores),
    }
    for name in metric_names:
        values = [float(row[name]) for row in per_target if is_number(row[name])]
        summary[f"mean_{name}"] = sum(values) / len(values) if values else math.nan
    return summary


def auroc(labels: list[int], scores: list[float]) -> float:
    positives = sum(labels)
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        return math.nan

    ranked = sorted(enumerate(scores), key=lambda item: item[1])
    ranks = [0.0] * len(scores)
    i = 0
    while i < len(ranked):
        j = i + 1
        while j < len(ranked) and ranked[j][1] == ranked[i][1]:
            j += 1
        average_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[ranked[k][0]] = average_rank
        i = j

    pos_rank_sum = sum(rank for rank, label in zip(ranks, labels) if label)
    return (pos_rank_sum - positives * (positives + 1) / 2.0) / (positives * negatives)


def average_precision(labels: list[int], scores: list[float]) -> float:
    positives = sum(labels)
    if positives == 0:
        return math.nan
    ordered = sorted(zip(scores, labels), key=lambda item: item[0], reverse=True)
    hits = 0
    precision_sum = 0.0
    for rank, (_, label) in enumerate(ordered, start=1):
        if label:
            hits += 1
            precision_sum += hits / rank
    return precision_sum / positives


def recall_at_k(labels: list[int], scores: list[float], k: int) -> float:
    positives = sum(labels)
    if positives == 0:
        return math.nan
    ordered = sorted(zip(scores, labels), key=lambda item: item[0], reverse=True)
    hits = sum(label for _, label in ordered[: max(0, k)])
    return hits / positives


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def print_summary(summary: dict[str, object], out_dir: Path) -> None:
    print(f"Targets: {summary['n_targets']}")
    print(f"Residues: {summary['n_residues']} ({summary['n_binding_residues']} binding)")
    print(f"Mean AUROC: {float(summary['mean_auroc']):.4f}")
    print(f"Mean AP: {float(summary['mean_average_precision']):.4f}")
    print(f"Mean F1 @ threshold: {float(summary['mean_f1']):.4f}")
    print(f"Mean recall @ threshold: {float(summary['mean_recall']):.4f}")
    print(f"Mean top-(N+2) recall: {float(summary['mean_top_nplus2_recall']):.4f}")
    print(f"Wrote evaluation tables: {out_dir}")


def residue_key(chain_id: object, residue_number: object) -> tuple[str, str]:
    return str(chain_id).strip(), str(residue_number).strip()


def split_codes(value: str) -> set[str]:
    return {item.strip() for item in value.split(",") if item.strip()}


def divide(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator else 0.0


def is_number(value: object) -> bool:
    try:
        return not math.isnan(float(value))
    except (TypeError, ValueError):
        return False


def resolve_path(value: str | Path, base: Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    return (base / path).resolve()


if __name__ == "__main__":
    raise SystemExit(main())
