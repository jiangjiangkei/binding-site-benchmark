#!/usr/bin/env python3
"""Run fpocket over a benchmark manifest and normalize residue scores."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO_ROOT / "data/manifests/p2rank_coach420_mlig_smoke20.csv"
DEFAULT_OUT_DIR = REPO_ROOT / "results/fpocket/p2rank_coach420_mlig_smoke20"
INFO_POCKET_RE = re.compile(r"^Pocket\s+(\d+)\s*:")
INFO_VALUE_RE = re.compile(r"^\s*([^:]+?)\s*:\s*([-+0-9.eE]+)")
HEADER_POCKET_RE = re.compile(r"Information about the pocket\s+(\d+)\s*:")
POCKET_FILE_RE = re.compile(r"pocket(\d+)_atm\.pdb$")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Batch-run fpocket prediction and write benchmark residue score TSV files.",
    )
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST), help="CSV manifest with target_id and structure_path.")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR), help="Directory for fpocket outputs.")
    parser.add_argument("--fpocket-bin", default="fpocket", help="fpocket executable.")
    parser.add_argument(
        "--probability-source",
        choices=("score", "druggability"),
        default="score",
        help="Pocket metric projected to residue probability. score is normalized per target.",
    )
    parser.add_argument("--threshold", type=float, default=0.5, help="Probability threshold recorded in is_binding.")
    parser.add_argument("--limit", type=int, help="Run at most this many manifest rows after filtering.")
    parser.add_argument("--target-id", action="append", default=[], help="Only run this target_id. Can be repeated.")
    parser.add_argument("--timeout", type=float, default=0.0, help="Per-target fpocket timeout in seconds. 0 disables timeout.")
    parser.add_argument("--resume", action="store_true", help="Skip targets whose fpocket raw output and score TSV exist.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running fpocket.")
    parser.add_argument("--fail-fast", action="store_true", help="Stop after the first failed target.")
    parser.add_argument("--extra-arg", action="append", default=[], help="Extra argument passed through to fpocket.")
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

    if not args.dry_run and not executable_exists(args.fpocket_bin):
        print(f"fpocket executable not found on PATH: {args.fpocket_bin}", file=sys.stderr)
        return 127

    (out_dir / "raw").mkdir(parents=True, exist_ok=True)
    (out_dir / "scores").mkdir(parents=True, exist_ok=True)
    (out_dir / "logs").mkdir(parents=True, exist_ok=True)

    summary_rows: list[dict[str, str]] = []
    exit_code = 0
    for row in rows:
        result = run_one(args, row, out_dir)
        summary_rows.append(result)
        print(f"{result['target_id']}: {result['status']} ({result['n_residues']} residues)")
        if result["status"] == "failed":
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
    raw_target_dir = out_dir / "raw" / target_id
    logs_dir = out_dir / "logs"
    scores_tsv = out_dir / "scores" / f"{target_id}.scores.tsv"
    work_pdb = raw_target_dir / f"{target_id}.pdb"
    fpocket_out_dir = raw_target_dir / f"{target_id}_out"
    info_path = fpocket_out_dir / f"{target_id}_info.txt"
    stdout_path = logs_dir / f"{target_id}.stdout.txt"
    stderr_path = logs_dir / f"{target_id}.stderr.txt"
    command = build_command(args, work_pdb)
    command_text = shlex.join(command)
    returncode = 0
    runtime = 0.0
    n_residues = 0

    if args.resume and info_path.exists() and scores_tsv.exists():
        n_residues = count_score_rows(scores_tsv)
        return summary(
            row,
            input_pdb,
            fpocket_out_dir,
            info_path,
            scores_tsv,
            n_residues,
            stdout_path,
            stderr_path,
            command_text,
            returncode,
            runtime,
            "skipped",
        )

    raw_target_dir.mkdir(parents=True, exist_ok=True)
    if not args.dry_run:
        shutil.copy2(input_pdb, work_pdb)

    if args.dry_run:
        print(command_text)
        status = "planned"
    else:
        started = time.monotonic()
        try:
            completed = subprocess.run(
                command,
                cwd=raw_target_dir,
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

        if returncode or not info_path.exists():
            status = "failed"
        else:
            try:
                n_residues = normalize_residue_scores(
                    input_pdb,
                    info_path,
                    fpocket_out_dir / "pockets",
                    scores_tsv,
                    row.get("chain_id", "").strip(),
                    args.threshold,
                    args.probability_source,
                )
                status = "ok"
            except ValueError as exc:
                print(f"{target_id}: {exc}", file=sys.stderr)
                status = "failed"

    return summary(
        row,
        input_pdb,
        fpocket_out_dir,
        info_path,
        scores_tsv,
        n_residues,
        stdout_path,
        stderr_path,
        command_text,
        returncode,
        runtime,
        status,
    )


def build_command(args: argparse.Namespace, work_pdb: Path) -> list[str]:
    command = [args.fpocket_bin, "-f", f"./{work_pdb.name}"]
    command.extend(args.extra_arg)
    return command


def normalize_residue_scores(
    input_pdb: Path,
    info_path: Path,
    pockets_dir: Path,
    scores_tsv: Path,
    chain_filter: str,
    threshold: float,
    probability_source: str,
) -> int:
    residue_rows = load_protein_residues(input_pdb, chain_filter)
    if chain_filter and not residue_rows:
        raise ValueError(f"No protein residues found for chain {chain_filter} in {input_pdb}")

    pocket_metrics = load_pocket_metrics(info_path)
    residue_scores = {residue_key(row["chain_id"], row["residue_number"]): [] for row in residue_rows}

    for pocket_path in sorted(pockets_dir.glob("pocket*_atm.pdb"), key=pocket_sort_key):
        pocket_number = pocket_number_from_file(pocket_path)
        if pocket_number is None:
            continue
        metrics = pocket_metrics.get(pocket_number, {})
        raw_score = metrics.get("score", 0.0)
        druggability = metrics.get("druggability_score", 0.0)
        for residue in residues_from_pocket_file(pocket_path, chain_filter):
            if residue in residue_scores:
                residue_scores[residue].append(
                    {
                        "pocket": pocket_number,
                        "fpocket_score": raw_score,
                        "druggability_score": druggability,
                    }
                )

    max_score = max((metrics.get("score", 0.0) for metrics in pocket_metrics.values()), default=0.0)
    rows = []
    for residue in residue_rows:
        key = residue_key(residue["chain_id"], residue["residue_number"])
        assignments = residue_scores[key]
        best = max(assignments, key=lambda item: item["fpocket_score"], default=None)
        if best is None:
            probability = 0.0
            fpocket_score = 0.0
            druggability_score = 0.0
            pocket = "0"
        else:
            fpocket_score = float(best["fpocket_score"])
            druggability_score = float(best["druggability_score"])
            pocket = str(best["pocket"])
            if probability_source == "druggability":
                probability = max(0.0, min(1.0, druggability_score))
            else:
                probability = fpocket_score / max_score if max_score > 0 else 0.0

        rows.append(
            {
                "residue_id": f"{residue['chain_id']}_{residue['residue_number']}",
                "chain_id": residue["chain_id"],
                "residue_number": residue["residue_number"],
                "probability": f"{probability:.6f}",
                "is_binding": "1" if probability >= threshold else "0",
                "fpocket_score": f"{fpocket_score:.6f}",
                "druggability_score": f"{druggability_score:.6f}",
                "pocket": pocket,
            }
        )

    scores_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "residue_id",
        "chain_id",
        "residue_number",
        "probability",
        "is_binding",
        "fpocket_score",
        "druggability_score",
        "pocket",
    ]
    with scores_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def load_protein_residues(path: Path, chain_filter: str) -> list[dict[str, str]]:
    seen = set()
    rows = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line[:6].strip() != "ATOM":
                continue
            atom = parse_pdb_atom(line)
            if atom is None:
                continue
            if chain_filter and atom["chain_id"] != chain_filter:
                continue
            key = residue_key(atom["chain_id"], atom["residue_number"])
            if key in seen:
                continue
            seen.add(key)
            rows.append({"chain_id": key[0], "residue_number": key[1]})
    return rows


def load_pocket_metrics(path: Path) -> dict[int, dict[str, float]]:
    if not path.exists():
        raise ValueError(f"Missing fpocket info file: {path}")

    metrics_by_pocket: dict[int, dict[str, float]] = {}
    current_pocket: int | None = None
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            pocket_match = INFO_POCKET_RE.match(line.strip())
            if pocket_match:
                current_pocket = int(pocket_match.group(1))
                metrics_by_pocket[current_pocket] = {}
                continue
            if current_pocket is None:
                continue
            value_match = INFO_VALUE_RE.match(line)
            if not value_match:
                continue
            key = normalize_metric_key(value_match.group(1))
            if key in {"score", "druggability_score"}:
                metrics_by_pocket[current_pocket][key] = float(value_match.group(2))
    return metrics_by_pocket


def residues_from_pocket_file(path: Path, chain_filter: str) -> set[tuple[str, str]]:
    residues = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line[:6].strip() != "ATOM":
                continue
            atom = parse_pdb_atom(line)
            if atom is None:
                continue
            if chain_filter and atom["chain_id"] != chain_filter:
                continue
            residues.add(residue_key(atom["chain_id"], atom["residue_number"]))
    return residues


def pocket_number_from_file(path: Path) -> int | None:
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = HEADER_POCKET_RE.search(line)
            if match:
                return int(match.group(1))

    match = POCKET_FILE_RE.search(path.name)
    if not match:
        return None
    pocket_index = int(match.group(1))
    return pocket_index + 1


def pocket_sort_key(path: Path) -> int:
    match = POCKET_FILE_RE.search(path.name)
    return int(match.group(1)) if match else 10**9


def parse_pdb_atom(line: str) -> dict[str, str] | None:
    if len(line) < 27:
        return None
    residue_number = line[22:26].strip()
    insertion_code = line[26].strip()
    if insertion_code:
        residue_number = f"{residue_number}{insertion_code}"
    return {
        "chain_id": line[21].strip(),
        "residue_number": residue_number,
    }


def normalize_metric_key(value: str) -> str:
    key = value.strip().lower()
    key = re.sub(r"[^a-z0-9]+", "_", key).strip("_")
    if key == "druggability_score" or key == "drug_score":
        return "druggability_score"
    return key


def summary(
    row: dict[str, str],
    input_pdb: Path,
    fpocket_out_dir: Path,
    info_path: Path,
    scores_tsv: Path,
    n_residues: int,
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
        "fpocket_out_dir": str(fpocket_out_dir),
        "fpocket_info": str(info_path),
        "scores_tsv": str(scores_tsv),
        "n_residues": str(n_residues),
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
        "fpocket_out_dir",
        "fpocket_info",
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


def count_score_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open(newline="", encoding="utf-8") as handle:
        return max(0, sum(1 for _ in handle) - 1)


def residue_key(chain_id: object, residue_number: object) -> tuple[str, str]:
    return str(chain_id).strip(), str(residue_number).strip()


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
