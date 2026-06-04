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
COLORS = ["#2563eb", "#dc2626", "#059669", "#7c3aed", "#ea580c", "#0891b2", "#9333ea", "#0f766e"]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot pooled residue-level ROC and PR curves.")
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    parser.add_argument("--scores-dir", default=str(DEFAULT_SCORES_DIR))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--score-suffix", default=".scores.tsv")
    parser.add_argument("--distance-cutoff", type=float, default=4.0)
    parser.add_argument("--method-name", default="ProtCross")
    parser.add_argument(
        "--method",
        action="append",
        default=[],
        help="Method curve as NAME=SCORES_DIR. Can be repeated. Uses --score-suffix for score files.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    manifest_path = ev.resolve_path(args.manifest, Path.cwd())
    scores_dir = ev.resolve_path(args.scores_dir, Path.cwd())
    out_dir = ev.resolve_path(args.out_dir, Path.cwd())
    out_dir.mkdir(parents=True, exist_ok=True)

    curves = []
    for method_name, method_scores_dir in parse_method_specs(args, scores_dir):
        labels, scores = collect_labels_and_scores(
            manifest_path,
            method_scores_dir,
            args.score_suffix,
            args.distance_cutoff,
        )
        curves.append(
            {
                "method_name": method_name,
                "scores_dir": str(method_scores_dir),
                "labels": labels,
                "scores": scores,
                "roc_points": roc_curve(labels, scores),
                "pr_points": precision_recall_curve(labels, scores),
                "auroc": ev.auroc(labels, scores),
                "average_precision": ev.average_precision(labels, scores),
                "positive_rate": sum(labels) / len(labels),
            }
        )

    if not curves:
        raise SystemExit("No method curves requested.")

    positive_rate = curves[0]["positive_rate"]

    roc_svg = out_dir / "roc_curve.svg"
    pr_svg = out_dir / "pr_curve.svg"
    roc_png = out_dir / "roc_curve.png"
    pr_png = out_dir / "pr_curve.png"
    metrics_csv = out_dir / "curve_metrics.csv"

    plot_roc(curves, roc_svg, roc_png)
    plot_pr(curves, positive_rate, pr_svg, pr_png)
    write_curve_metrics(metrics_csv, curves)

    print(f"Wrote ROC curve: {roc_svg}")
    print(f"Wrote PR curve: {pr_svg}")
    print(f"Wrote curve metrics: {metrics_csv}")
    for curve in curves:
        print(f"{curve['method_name']}: AUROC {curve['auroc']:.4f}; AP {curve['average_precision']:.4f}")
    return 0


def parse_method_specs(args: argparse.Namespace, legacy_scores_dir: Path) -> list[tuple[str, Path]]:
    if not args.method:
        return [(args.method_name, legacy_scores_dir)]

    specs = []
    for value in args.method:
        if "=" not in value:
            raise SystemExit(f"Invalid --method value, expected NAME=SCORES_DIR: {value}")
        name, scores_dir = value.split("=", 1)
        name = name.strip()
        scores_dir = scores_dir.strip()
        if not name or not scores_dir:
            raise SystemExit(f"Invalid --method value, expected NAME=SCORES_DIR: {value}")
        specs.append((name, ev.resolve_path(scores_dir, Path.cwd())))
    return specs


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
    curves: list[dict[str, object]],
    svg_path: Path,
    png_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 5.2), dpi=160)
    style_axes(ax)
    ax.plot([0, 1], [0, 1], color="#9ca3af", linestyle="--", linewidth=1.2, label="Random")
    for index, curve in enumerate(curves):
        x, y = zip(*curve["roc_points"])
        color = COLORS[index % len(COLORS)]
        ax.plot(
            x,
            y,
            color=color,
            linewidth=2.2,
            label=f"{curve['method_name']} AUROC = {curve['auroc']:.3f}",
        )
    ax.set_title("Residue-Level ROC")
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
    curves: list[dict[str, object]],
    positive_rate: float,
    svg_path: Path,
    png_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 5.2), dpi=160)
    style_axes(ax)
    ax.axhline(
        positive_rate,
        color="#9ca3af",
        linestyle="--",
        linewidth=1.2,
        label=f"Baseline = {positive_rate:.3f}",
    )
    for index, curve in enumerate(curves):
        x, y = zip(*curve["pr_points"])
        color = COLORS[index % len(COLORS)]
        ax.plot(
            x,
            y,
            color=color,
            linewidth=2.2,
            label=f"{curve['method_name']} AP = {curve['average_precision']:.3f}",
        )
    ax.set_title("Residue-Level PR")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(loc="upper right", frameon=False)
    fig.tight_layout()
    fig.savefig(svg_path)
    fig.savefig(png_path)
    plt.close(fig)


def write_curve_metrics(path: Path, curves: list[dict[str, object]]) -> None:
    fieldnames = [
        "method_name",
        "scores_dir",
        "n_residues",
        "n_binding_residues",
        "positive_rate",
        "auroc",
        "average_precision",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for curve in curves:
            labels = curve["labels"]
            writer.writerow(
                {
                    "method_name": curve["method_name"],
                    "scores_dir": curve["scores_dir"],
                    "n_residues": len(labels),
                    "n_binding_residues": sum(labels),
                    "positive_rate": curve["positive_rate"],
                    "auroc": curve["auroc"],
                    "average_precision": curve["average_precision"],
                }
            )


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
