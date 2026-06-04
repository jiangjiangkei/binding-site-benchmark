#!/usr/bin/env python3
"""Run P2Rank over a benchmark manifest and normalize residue scores."""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO_ROOT / "data/manifests/p2rank_coach420_mlig_smoke20.csv"
DEFAULT_OUT_DIR = REPO_ROOT / "results/p2rank/p2rank_coach420_mlig_smoke20"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Batch-run P2Rank prediction and write benchmark residue score TSV files.",
    )
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST), help="CSV manifest with target_id and structure_path.")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR), help="Directory for P2Rank outputs.")
    parser.add_argument("--p2rank-bin", default="prank", help="P2Rank prank executable.")
    parser.add_argument("--config", help="P2Rank config passed with -c, for example alphafold.")
    parser.add_argument("--model", help="P2Rank model passed with -m.")
    parser.add_argument("--threads", type=int, default=1, help="Threads passed to P2Rank.")
    parser.add_argument("--visualizations", choices=("0", "1"), default="0", help="Whether P2Rank writes visualizations.")
    parser.add_argument("--threshold", type=float, default=0.5, help="Probability threshold recorded in is_binding.")
    parser.add_argument("--limit", type=int, help="Run at most this many manifest rows after filtering.")
    parser.add_argument("--target-id", action="append", default=[], help="Only run this target_id. Can be repeated.")
    parser.add_argument("--timeout", type=float, default=0.0, help="Whole P2Rank run timeout in seconds. 0 disables timeout.")
    parser.add_argument("--resume", action="store_true", help="Skip P2Rank prediction for targets with raw residue CSVs.")
    parser.add_argument("--dry-run", action="store_true", help="Print the P2Rank command without running it.")
    parser.add_argument("--extra-arg", action="append", default=[], help="Extra argument passed through to P2Rank.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest_path = resolve_path(args.manifest, Path.cwd())
    out_dir = resolve_path(args.out_dir, Path.cwd())
    rows = load_manifest(manifest_path, set(args.target_id))
    if args.limit is not None:
        rows = rows[: args.limit]

    if not rows:
        print("No manifest rows selected.", file=sys.stderr)
        return 2

    if not args.dry_run and not executable_exists(args.p2rank_bin):
        print(f"P2Rank executable not found on PATH: {args.p2rank_bin}", file=sys.stderr)
        return 127

    raw_dir = out_dir / "raw"
    scores_dir = out_dir / "scores"
    logs_dir = out_dir / "logs"
    raw_dir.mkdir(parents=True, exist_ok=True)
    scores_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    rows_to_predict = rows
    if args.resume:
        rows_to_predict = [row for row in rows if not p2rank_residues_path(row, raw_dir).exists()]

    command: list[str] = []
    command_text = ""
    returncode = 0
    runtime = 0.0
    stdout_path = logs_dir / "p2rank.stdout.txt"
    stderr_path = logs_dir / "p2rank.stderr.txt"
    dataset_path = out_dir / "p2rank_input.ds"

    if rows_to_predict:
        write_prediction_dataset(dataset_path, rows_to_predict)
        command = build_command(args, dataset_path, raw_dir)
        command_text = shlex.join(command)

        if args.dry_run:
            print(command_text)
        else:
            started = time.monotonic()
            try:
                completed = subprocess.run(
                    command,
                    cwd=REPO_ROOT,
                    capture_output=True,
                    text=True,
                    timeout=args.timeout or None,
                    check=False,
                )
                runtime = time.monotonic() - started
                stdout_path.write_text(completed.stdout, encoding="utf-8")
                stderr_path.write_text(completed.stderr, encoding="utf-8")
                returncode = completed.returncode
            except subprocess.TimeoutExpired as exc:
                runtime = time.monotonic() - started
                stdout_path.write_text(exc.stdout or "", encoding="utf-8")
                stderr_path.write_text((exc.stderr or "") + f"\nTimed out after {args.timeout} seconds.\n", encoding="utf-8")
                returncode = 124
    elif args.resume:
        print("All selected targets have P2Rank raw residue CSVs; normalizing scores.")

    summary_rows = summarize_outputs(
        rows,
        raw_dir,
        scores_dir,
        stdout_path,
        stderr_path,
        command_text,
        returncode,
        runtime,
        args.threshold,
        args.dry_run,
        {row["target_id"] for row in rows_to_predict},
        args.resume,
    )
    write_summary(out_dir / "run_summary.csv", summary_rows)

    failed = [row for row in summary_rows if row["status"] == "failed"]
    for row in summary_rows:
        print(f"{row['target_id']}: {row['status']} ({row['n_residues']} residues)")
    print(f"Wrote summary: {out_dir / 'run_summary.csv'}")
    return 1 if returncode or failed else 0


def load_manifest(path: Path, selected_target_ids: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"target_id", "structure_path"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Manifest missing required columns: {', '.join(sorted(missing))}")

        rows = []
        for row in reader:
            if selected_target_ids and row["target_id"] not in selected_target_ids:
                continue
            rows.append(row)
        return rows


def build_command(args: argparse.Namespace, dataset_path: Path, raw_dir: Path) -> list[str]:
    command = [
        args.p2rank_bin,
        "predict",
        str(dataset_path),
        "-o",
        str(raw_dir),
        "-threads",
        str(args.threads),
        "-visualizations",
        args.visualizations,
    ]
    if args.config:
        command.extend(["-c", args.config])
    if args.model:
        command.extend(["-m", args.model])
    command.extend(args.extra_arg)
    return command


def write_prediction_dataset(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(str(resolve_path(row["structure_path"], REPO_ROOT)))
            handle.write("\n")


def summarize_outputs(
    rows: list[dict[str, str]],
    raw_dir: Path,
    scores_dir: Path,
    stdout_path: Path,
    stderr_path: Path,
    command_text: str,
    returncode: int,
    runtime: float,
    threshold: float,
    dry_run: bool,
    predicted_target_ids: set[str],
    resume: bool,
) -> list[dict[str, str]]:
    summary_rows = []
    for row in rows:
        target_id = safe_name(row["target_id"])
        input_pdb = resolve_path(row["structure_path"], REPO_ROOT)
        p2rank_csv = p2rank_residues_path(row, raw_dir)
        scores_tsv = scores_dir / f"{target_id}.scores.tsv"
        n_residues = 0

        if dry_run:
            status = "planned"
        elif returncode:
            status = "failed"
        elif not p2rank_csv.exists():
            status = "failed"
        else:
            try:
                score_existed = scores_tsv.exists()
                if not (resume and score_existed and row["target_id"] not in predicted_target_ids):
                    n_residues = normalize_residue_scores(
                        p2rank_csv,
                        scores_tsv,
                        row.get("chain_id", "").strip(),
                        threshold,
                    )
                else:
                    n_residues = count_score_rows(scores_tsv)
                status = "skipped" if resume and score_existed and row["target_id"] not in predicted_target_ids else "ok"
            except ValueError as exc:
                print(f"{target_id}: {exc}", file=sys.stderr)
                status = "failed"

        summary_rows.append(
            {
                "target_id": row.get("target_id", ""),
                "pdb_id": row.get("pdb_id", ""),
                "chain_id": row.get("chain_id", ""),
                "input_pdb": str(input_pdb),
                "p2rank_residues_csv": str(p2rank_csv),
                "scores_tsv": str(scores_tsv),
                "n_residues": str(n_residues),
                "status": status,
                "returncode": str(returncode),
                "runtime_sec": f"{runtime:.3f}",
                "stdout_path": str(stdout_path),
                "stderr_path": str(stderr_path),
                "command": command_text,
            }
        )
    return summary_rows


def normalize_residue_scores(p2rank_csv: Path, scores_tsv: Path, chain_filter: str, threshold: float) -> int:
    rows = []
    with p2rank_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, skipinitialspace=True)
        fieldnames = {field.strip(): field for field in reader.fieldnames or []}
        required = {"chain", "residue_label", "score", "probability", "pocket"}
        missing = required.difference(fieldnames)
        if missing:
            raise ValueError(f"P2Rank residue CSV missing columns: {p2rank_csv}: {', '.join(sorted(missing))}")

        for row in reader:
            chain_id = row[fieldnames["chain"]].strip()
            if chain_filter and chain_id != chain_filter:
                continue
            residue_number = row[fieldnames["residue_label"]].strip()
            probability = float(row[fieldnames["probability"]])
            score = float(row[fieldnames["score"]])
            pocket = row[fieldnames["pocket"]].strip()
            rows.append(
                {
                    "residue_id": f"{chain_id}_{residue_number}",
                    "chain_id": chain_id,
                    "residue_number": residue_number,
                    "probability": f"{probability:.6f}",
                    "is_binding": "1" if probability >= threshold else "0",
                    "p2rank_score": f"{score:.6f}",
                    "pocket": pocket,
                }
            )

    if chain_filter and not rows:
        raise ValueError(f"No P2Rank residues found for chain {chain_filter} in {p2rank_csv}")

    scores_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames_out = ["residue_id", "chain_id", "residue_number", "probability", "is_binding", "p2rank_score", "pocket"]
    with scores_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames_out, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def count_score_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open(newline="", encoding="utf-8") as handle:
        return max(0, sum(1 for _ in handle) - 1)


def p2rank_residues_path(row: dict[str, str], raw_dir: Path) -> Path:
    input_pdb = resolve_path(row["structure_path"], REPO_ROOT)
    return raw_dir / f"{input_pdb.name}_residues.csv"


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "target_id",
        "pdb_id",
        "chain_id",
        "input_pdb",
        "p2rank_residues_csv",
        "scores_tsv",
        "n_residues",
        "status",
        "returncode",
        "runtime_sec",
        "stdout_path",
        "stderr_path",
        "command",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def executable_exists(value: str) -> bool:
    return shutil.which(value) is not None or Path(value).expanduser().exists()


def resolve_path(value: str | os.PathLike[str], base: Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    return (base / path).resolve()


def safe_name(value: str) -> str:
    return "".join(char if char.isalnum() or char in ("-", "_", ".") else "_" for char in value)


if __name__ == "__main__":
    raise SystemExit(main())
