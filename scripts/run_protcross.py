#!/usr/bin/env python3
"""Run ProtCross over a benchmark manifest.

The script is intentionally a thin adapter around the official ProtCross CLI.
It keeps this repository's benchmark bookkeeping separate from ProtCross model
code and writes one residue score table per target.
"""

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
DEFAULT_OUT_DIR = REPO_ROOT / "results/protcross/p2rank_coach420_mlig_smoke20"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Batch-run ProtCross prediction for structures listed in a benchmark manifest.",
    )
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST), help="CSV manifest with target_id and structure_path.")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR), help="Directory for ProtCross outputs.")
    parser.add_argument("--protcross-bin", default="protcross", help="ProtCross CLI executable.")
    parser.add_argument("--assets-dir", help="Directory containing ProtCross checkpoint, PCA, and ESM-C weights.")
    parser.add_argument("--checkpoint", help="Explicit ProtCross checkpoint path.")
    parser.add_argument("--esm-weights", help="Explicit ESM-C weights path.")
    parser.add_argument("--pca", help="Explicit ProtCross PCA reducer path.")
    parser.add_argument("--threshold", type=float, default=0.5, help="Residue probability threshold recorded by ProtCross.")
    parser.add_argument("--device", default="auto", help="ProtCross device: auto, cpu, cuda, or cuda:N.")
    parser.add_argument("--limit", type=int, help="Run at most this many manifest rows after filtering.")
    parser.add_argument("--target-id", action="append", default=[], help="Only run this target_id. Can be repeated.")
    parser.add_argument("--timeout", type=float, default=0.0, help="Per-target timeout in seconds. 0 disables timeout.")
    parser.add_argument("--resume", action="store_true", help="Skip targets whose output PDB and score TSV already exist.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running ProtCross.")
    parser.add_argument("--fail-fast", action="store_true", help="Stop after the first failed target.")
    parser.add_argument("--extra-arg", action="append", default=[], help="Extra argument passed through to ProtCross.")
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

    if not args.dry_run and shutil.which(args.protcross_bin) is None:
        print(f"ProtCross executable not found on PATH: {args.protcross_bin}", file=sys.stderr)
        return 127

    (out_dir / "pdb").mkdir(parents=True, exist_ok=True)
    (out_dir / "scores").mkdir(parents=True, exist_ok=True)
    (out_dir / "logs").mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, str]] = []
    exit_code = 0
    for row in rows:
        result = run_one(args, row, out_dir)
        summary_rows.append(result)
        status = result["status"]
        print(f"{result['target_id']}: {status} ({result['runtime_sec']}s)")
        if status == "failed":
            exit_code = 1
            if args.fail_fast:
                break

    write_summary(out_dir / "run_summary.csv", summary_rows)
    print(f"Wrote summary: {out_dir / 'run_summary.csv'}")
    return exit_code


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


def run_one(args: argparse.Namespace, row: dict[str, str], out_dir: Path) -> dict[str, str]:
    target_id = safe_name(row["target_id"])
    input_pdb = resolve_path(row["structure_path"], REPO_ROOT)
    output_pdb = out_dir / "pdb" / f"{target_id}.protcross.pdb"
    scores_tsv = out_dir / "scores" / f"{target_id}.scores.tsv"
    stdout_path = out_dir / "logs" / f"{target_id}.stdout.txt"
    stderr_path = out_dir / "logs" / f"{target_id}.stderr.txt"

    command = build_command(args, row, input_pdb, output_pdb, scores_tsv)
    command_text = shlex.join(command)
    started = time.monotonic()

    if args.resume and output_pdb.exists() and scores_tsv.exists():
        return summary(row, input_pdb, output_pdb, scores_tsv, stdout_path, stderr_path, command_text, 0, 0.0, "skipped")

    if args.dry_run:
        print(command_text)
        return summary(row, input_pdb, output_pdb, scores_tsv, stdout_path, stderr_path, command_text, 0, 0.0, "planned")

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
        status = "ok" if completed.returncode == 0 else "failed"
        return summary(
            row,
            input_pdb,
            output_pdb,
            scores_tsv,
            stdout_path,
            stderr_path,
            command_text,
            completed.returncode,
            runtime,
            status,
        )
    except subprocess.TimeoutExpired as exc:
        runtime = time.monotonic() - started
        stdout_path.write_text(exc.stdout or "", encoding="utf-8")
        stderr_path.write_text((exc.stderr or "") + f"\nTimed out after {args.timeout} seconds.\n", encoding="utf-8")
        return summary(row, input_pdb, output_pdb, scores_tsv, stdout_path, stderr_path, command_text, 124, runtime, "failed")


def build_command(
    args: argparse.Namespace,
    row: dict[str, str],
    input_pdb: Path,
    output_pdb: Path,
    scores_tsv: Path,
) -> list[str]:
    executable = args.protcross_bin
    command = [executable]
    if Path(executable).name != "protcross-predict":
        command.append("predict")

    command.extend(
        [
            str(input_pdb),
            "--output",
            str(output_pdb),
            "--scores-tsv",
            str(scores_tsv),
            "--threshold",
            str(args.threshold),
            "--device",
            args.device,
            "--quiet",
        ]
    )

    chain_id = row.get("chain_id", "").strip()
    if chain_id:
        command.extend(["--chain", chain_id])

    optional_paths = [
        ("--assets-dir", args.assets_dir),
        ("--checkpoint", args.checkpoint),
        ("--esm-weights", args.esm_weights),
        ("--pca", args.pca),
    ]
    for flag, value in optional_paths:
        if value:
            command.extend([flag, str(resolve_path(value, Path.cwd()))])

    command.extend(args.extra_arg)
    return command


def summary(
    row: dict[str, str],
    input_pdb: Path,
    output_pdb: Path,
    scores_tsv: Path,
    stdout_path: Path,
    stderr_path: Path,
    command: str,
    returncode: int,
    runtime: float,
    status: str,
) -> dict[str, str]:
    return {
        "target_id": row.get("target_id", ""),
        "pdb_id": row.get("pdb_id", ""),
        "chain_id": row.get("chain_id", ""),
        "input_pdb": str(input_pdb),
        "output_pdb": str(output_pdb),
        "scores_tsv": str(scores_tsv),
        "status": status,
        "returncode": str(returncode),
        "runtime_sec": f"{runtime:.3f}",
        "stdout_path": str(stdout_path),
        "stderr_path": str(stderr_path),
        "command": command,
    }


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "target_id",
        "pdb_id",
        "chain_id",
        "input_pdb",
        "output_pdb",
        "scores_tsv",
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


def resolve_path(value: str | os.PathLike[str], base: Path) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    return (base / path).resolve()


def safe_name(value: str) -> str:
    return "".join(char if char.isalnum() or char in ("-", "_", ".") else "_" for char in value)


if __name__ == "__main__":
    raise SystemExit(main())
