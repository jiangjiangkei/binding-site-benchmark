#!/usr/bin/env python3
"""Plot pooled residue-level ROC and PR curves."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt

import evaluate_residue_scores as ev


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO_ROOT / "data/manifests/p2rank_coach420_mlig_smoke20.csv"
DEFAULT_SCORES_DIR = REPO_ROOT / "results/protcross/p2rank_coach420_mlig_smoke20/scores"
DEFAULT_OUT_DIR = REPO_ROOT / "results/protcross/p2rank_coach420_mlig_smoke20/evaluation_residue"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot pooled residue-level ROC and PR curves.")
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--scores-dir", default=str(DEFAULT_SCORES_DIR))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--score-suffix", default=".scores.tsv")
    parser.add_argument("--distance-cutoff", type=float, default=4.0)
    parser.add_argument("--method-name", default="ProtCross")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    manifest_path = ev.resolve_path(args.manifest, Path.cwd())
    scores_dir = ev.resolve_path(args.scores_dir, Path.cwd())
    out_dir = ev.resolve_path(args.out_dir, Path.cwd())
    out_dir.mkdir(parents=True, exist_ok=True)

    labels, scores = collect_labels_and_scores(manifest_path, scores_dir, args.score_suffix, args.distance_cutoff)
    roc_points = roc_curve(labels, scores)
    pr_points = precision_recall_curve(labels, scores)
    auroc = ev.auroc(labels, scores)
    average_precision = ev.average_precision(labels, scores)
    positive_rate = sum(labels) / len(labels)

    roc_svg = out_dir / "roc_curve.svg"
    pr_svg = out_dir / "pr_curve.svg"
    roc_png = out_dir / "roc_curve.png"
    pr_png = out_dir / "pr_curve.png"

    plot_roc(roc_points, auroc, args.method_name, roc_svg, roc_png)
    plot_pr(pr_points, average_precision, positive_rate, args.method_name, pr_svg, pr_png)

    print(f"Wrote ROC curve: {roc_svg}")
    print(f"Wrote PR curve: {pr_svg}")
    print(f"AUROC: {auroc:.4f}")
    print(f"Average precision: {average_precision:.4f}")
    return 0


def collect_labels_and_scores(
    manifest_path: Path,
    scores_dir: Path,
    score_suffix: str,
    distance_cutoff: float,
) -> tuple[list[int], list[float]]:
    labels: list[int] = []
    scores: list[float] = []
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            target_id = row["target_id"]
            score_records = ev.load_scores(scores_dir / f"{target_id}{score_suffix}")
            truth = ev.binding_residues_from_pdb(
                ev.resolve_path(row["structure_path"], REPO_ROOT),
                ev.split_codes(row["ligand_codes"]),
                row.get("chain_id", "").strip(),
                distance_cutoff,
            )
            for record in score_records:
                labels.append(1 if ev.residue_key(record["chain_id"], record["residue_number"]) in truth else 0)
                scores.append(float(record["probability"]))
    return labels, scores


def roc_curve(labels: list[int], scores: list[float]) -> list[tuple[float, float]]:
    positives = sum(labels)
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        raise SystemExit("ROC requires at least one positive and one negative residue.")

    ordered = sorted(zip(scores, labels), key=lambda item: item[0], reverse=True)
    points = [(0.0, 0.0)]
    true_positive = 0
    false_positive = 0
    last_score = None
    for score, label in ordered:
        if last_score is not None and score != last_score:
            points.append((false_positive / negatives, true_positive / positives))
        if label:
            true_positive += 1
        else:
            false_positive += 1
        last_score = score
    points.append((false_positive / negatives, true_positive / positives))
    return points


def precision_recall_curve(labels: list[int], scores: list[float]) -> list[tuple[float, float]]:
    positives = sum(labels)
    if positives == 0:
        raise SystemExit("PR curve requires at least one positive residue.")

    ordered = sorted(zip(scores, labels), key=lambda item: item[0], reverse=True)
    points = [(0.0, 1.0)]
    true_positive = 0
    false_positive = 0
    last_score = None
    for score, label in ordered:
        if last_score is not None and score != last_score:
            points.append((true_positive / positives, precision(true_positive, false_positive)))
        if label:
            true_positive += 1
        else:
            false_positive += 1
        last_score = score
    points.append((true_positive / positives, precision(true_positive, false_positive)))
    return points


def precision(true_positive: int, false_positive: int) -> float:
    denominator = true_positive + false_positive
    return true_positive / denominator if denominator else 1.0


def plot_roc(
    points: list[tuple[float, float]],
    auroc: float,
    method_name: str,
    svg_path: Path,
    png_path: Path,
) -> None:
    x, y = zip(*points)
    fig, ax = plt.subplots(figsize=(6.2, 5.2), dpi=160)
    style_axes(ax)
    ax.plot([0, 1], [0, 1], color="#9ca3af", linestyle="--", linewidth=1.2, label="Random")
    ax.plot(x, y, color="#2563eb", linewidth=2.2, label=f"{method_name} AUROC = {auroc:.3f}")
    ax.set_title(f"{method_name} Residue-Level ROC")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(loc="lower right", frameon=False)
    fig.tight_layout()
    fig.savefig(svg_path)
    fig.savefig(png_path)
    plt.close(fig)


def plot_pr(
    points: list[tuple[float, float]],
    average_precision: float,
    positive_rate: float,
    method_name: str,
    svg_path: Path,
    png_path: Path,
) -> None:
    x, y = zip(*points)
    fig, ax = plt.subplots(figsize=(6.2, 5.2), dpi=160)
    style_axes(ax)
    ax.axhline(
        positive_rate,
        color="#9ca3af",
        linestyle="--",
        linewidth=1.2,
        label=f"Baseline = {positive_rate:.3f}",
    )
    ax.plot(x, y, color="#dc2626", linewidth=2.2, label=f"{method_name} AP = {average_precision:.3f}")
    ax.set_title(f"{method_name} Residue-Level PR")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    fig.savefig(svg_path)
    fig.savefig(png_path)
    plt.close(fig)


def style_axes(ax: plt.Axes) -> None:
    ax.set_facecolor("#ffffff")
    ax.grid(True, color="#e5e7eb", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#6b7280")
    ax.spines["bottom"].set_color("#6b7280")
    ax.tick_params(colors="#374151")
    ax.title.set_color("#111827")
    ax.xaxis.label.set_color("#111827")
    ax.yaxis.label.set_color("#111827")


if __name__ == "__main__":
    raise SystemExit(main())
