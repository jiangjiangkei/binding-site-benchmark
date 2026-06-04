#!/usr/bin/env python3
"""Run Kalasanty over a benchmark manifest and normalize residue scores."""

from __future__ import annotations

import argparse
import csv
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

import residue_score_utils as rsu


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO_ROOT / "data/manifests/p2rank_coach420_mlig_smoke20.csv"
DEFAULT_OUT_DIR = REPO_ROOT / "results/kalasanty/p2rank_coach420_mlig_smoke20"
POCKET_MOL2_RE = re.compile(r"pocket(\d+)\.mol2$")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Batch-run Kalasanty prediction and write benchmark residue score TSV files.",
    )
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST), help="CSV manifest with target_id and structure_path.")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR), help="Directory for Kalasanty outputs.")
    parser.add_argument("--kalasanty-dir", help="Directory containing Kalasanty scripts/predict.py and data/model_scpdb2017.hdf.")
    parser.add_argument("--predict-script", help="Explicit path to Kalasanty scripts/predict.py.")
    parser.add_argument("--model", help="Explicit Kalasanty .hdf model path.")
    parser.add_argument("--python-bin", default=sys.executable, help="Python executable in the Kalasanty environment.")
    parser.add_argument("--docker-image", help="Run Kalasanty through this Docker image instead of --python-bin.")
    parser.add_argument("--docker-platform", default="linux/amd64", help="Docker platform used with --docker-image.")
    parser.add_argument("--docker-sudo", action="store_true", help="Prefix Docker commands with sudo -n.")
    parser.add_argument("--docker-kalasanty-dir", default="/opt/kalasanty", help="Kalasanty directory inside --docker-image.")
    parser.add_argument("--input-format", default="pdb", help="Input format passed to Kalasanty/Open Babel.")
    parser.add_argument("--threshold", type=float, default=0.5, help="Residue probability threshold recorded in is_binding.")
    parser.add_argument("--coordinate-tolerance", type=float, default=0.25, help="Max Angstrom distance to map MOL2 atoms back to PDB residues.")
    parser.add_argument(
        "--pocket-score-mode",
        choices=("binary", "rank-decay"),
        default="binary",
        help="Kalasanty pocket membership scoring. binary uses 1.0 for all predicted pockets.",
    )
    parser.add_argument("--limit", type=int, help="Run at most this many manifest rows after filtering.")
    parser.add_argument("--target-id", action="append", default=[], help="Only run this target_id. Can be repeated.")
    parser.add_argument("--timeout", type=float, default=0.0, help="Per-target Kalasanty timeout in seconds. 0 disables timeout.")
    parser.add_argument("--resume", action="store_true", help="Skip targets whose raw output and score TSV already exist.")
    parser.add_argument("--normalize-only", action="store_true", help="Only parse existing Kalasanty raw outputs.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running Kalasanty.")
    parser.add_argument("--fail-fast", action="store_true", help="Stop after the first failed target.")
    parser.add_argument("--extra-arg", action="append", default=[], help="Extra argument passed through to Kalasanty predict.py.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest_path = rsu.resolve_path(args.manifest, Path.cwd())
    out_dir = rsu.resolve_path(args.out_dir, Path.cwd())
    rows = rsu.load_manifest(manifest_path, set(args.target_id))
    if args.limit is not None:
        rows = rows[: args.limit]

    if not rows:
        print("No manifest rows selected.", file=sys.stderr)
        return 2

    if not args.dry_run and not args.normalize_only:
        validate_runtime_args(args)

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


def validate_runtime_args(args: argparse.Namespace) -> None:
    if args.docker_image:
        if shutil.which("docker") is None:
            raise SystemExit("Docker executable not found for --docker-image run.")
        return

    predict_script = kalasanty_predict_script(args)
    if not predict_script.exists():
        raise SystemExit(f"Kalasanty predict.py not found: {predict_script}")
    if args.model and not rsu.resolve_path(args.model, Path.cwd()).exists():
        raise SystemExit(f"Kalasanty model not found: {args.model}")
    if not rsu.executable_exists(args.python_bin):
        raise SystemExit(f"Python executable not found: {args.python_bin}")


def run_one(args: argparse.Namespace, row: dict[str, str], out_dir: Path) -> dict[str, str]:
    target_id = rsu.safe_name(row["target_id"])
    input_pdb = rsu.resolve_path(row["structure_path"], REPO_ROOT)
    raw_target_dir = out_dir / "raw" / target_id
    logs_dir = out_dir / "logs"
    scores_tsv = out_dir / "scores" / f"{target_id}.scores.tsv"
    work_pdb = raw_target_dir / f"{target_id}.protein.pdb"
    prediction_root = raw_target_dir / "predicted_pockets"
    prediction_dir = prediction_root / target_id
    stdout_path = logs_dir / f"{target_id}.stdout.txt"
    stderr_path = logs_dir / f"{target_id}.stderr.txt"
    command = build_command(args, work_pdb, prediction_root)
    command_text = shlex.join(command)
    returncode = 0
    runtime = 0.0
    n_residues = 0

    if args.resume and pocket_mol2_files(prediction_dir) and scores_tsv.exists():
        n_residues = rsu.count_score_rows(scores_tsv)
        return summary(row, input_pdb, prediction_dir, scores_tsv, n_residues, stdout_path, stderr_path, command_text, 0, 0.0, "skipped")

    raw_target_dir.mkdir(parents=True, exist_ok=True)
    if not args.dry_run:
        rsu.write_protein_only_pdb(input_pdb, work_pdb, row.get("chain_id", "").strip())

    if args.dry_run:
        print(command_text)
        status = "planned"
    else:
        if not args.normalize_only:
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

        if returncode:
            status = "failed"
        else:
            try:
                n_residues = normalize_residue_scores(
                    input_pdb,
                    prediction_dir,
                    scores_tsv,
                    row.get("chain_id", "").strip(),
                    args.threshold,
                    args.coordinate_tolerance,
                    args.pocket_score_mode,
                )
                status = "ok"
            except ValueError as exc:
                print(f"{target_id}: {exc}", file=sys.stderr)
                status = "failed"

    return summary(row, input_pdb, prediction_dir, scores_tsv, n_residues, stdout_path, stderr_path, command_text, returncode, runtime, status)


def build_command(args: argparse.Namespace, work_pdb: Path, prediction_root: Path) -> list[str]:
    if args.docker_image:
        command = []
        if args.docker_sudo:
            command.extend(["sudo", "-n"])
        command.extend(
            [
                "docker",
                "run",
                "--rm",
                "--platform",
                args.docker_platform,
                "-v",
                f"{REPO_ROOT.parent}:{REPO_ROOT.parent}",
                "-v",
                "/tmp:/tmp",
                "-w",
                str(work_pdb.parent),
                args.docker_image,
                str(Path(args.docker_kalasanty_dir) / "scripts" / "predict.py"),
                "--input",
                str(work_pdb),
                "--output",
                str(prediction_root),
                "--format",
                args.input_format,
                "--dirname_pattern",
                r".*/([^/]+)/[^/]+$",
            ]
        )
        if args.model:
            command.extend(["--model", str(rsu.resolve_path(args.model, Path.cwd()))])
        command.extend(args.extra_arg)
        return command

    command = [
        args.python_bin,
        str(kalasanty_predict_script(args)),
        "--input",
        str(work_pdb),
        "--output",
        str(prediction_root),
        "--format",
        args.input_format,
        "--dirname_pattern",
        r".*/([^/]+)/[^/]+$",
    ]
    if args.model:
        command.extend(["--model", str(rsu.resolve_path(args.model, Path.cwd()))])
    command.extend(args.extra_arg)
    return command


def kalasanty_predict_script(args: argparse.Namespace) -> Path:
    if args.predict_script:
        return rsu.resolve_path(args.predict_script, Path.cwd())
    if not args.kalasanty_dir:
        return Path("scripts/predict.py")
    return rsu.resolve_path(args.kalasanty_dir, Path.cwd()) / "scripts" / "predict.py"


def normalize_residue_scores(
    input_pdb: Path,
    prediction_dir: Path,
    scores_tsv: Path,
    chain_filter: str,
    threshold: float,
    coordinate_tolerance: float,
    pocket_score_mode: str,
) -> int:
    residue_rows = rsu.load_protein_residues(input_pdb, chain_filter)
    protein_atoms = rsu.load_protein_atoms(input_pdb, chain_filter)
    if chain_filter and not residue_rows:
        raise ValueError(f"No protein residues found for chain {chain_filter} in {input_pdb}")

    residue_scores: dict[tuple[str, str], dict[str, str | float]] = {}
    pockets = pocket_mol2_files(prediction_dir)
    for rank_index, pocket_path in enumerate(pockets, start=1):
        coordinates = rsu.mol2_atom_coordinates(pocket_path)
        residues = rsu.residues_from_coordinates(coordinates, protein_atoms, coordinate_tolerance)
        probability = 1.0 if pocket_score_mode == "binary" else 1.0 / rank_index
        pocket = str(pocket_number(pocket_path))
        for residue in residues:
            current = residue_scores.get(residue)
            if current is None or probability > float(current["probability"]):
                residue_scores[residue] = {"probability": probability, "pocket": pocket}

    return rsu.write_residue_scores(scores_tsv, residue_rows, residue_scores, "kalasanty_score", threshold)


def pocket_mol2_files(prediction_dir: Path) -> list[Path]:
    if not prediction_dir.exists():
        return []
    return sorted(prediction_dir.glob("pocket*.mol2"), key=pocket_number)


def pocket_number(path: Path) -> int:
    match = POCKET_MOL2_RE.match(path.name)
    return int(match.group(1)) if match else 10**9


def summary(
    row: dict[str, str],
    input_pdb: Path,
    prediction_dir: Path,
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
        "prediction_dir": str(prediction_dir),
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
        "prediction_dir",
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


if __name__ == "__main__":
    raise SystemExit(main())
