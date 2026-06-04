#!/usr/bin/env python3
"""Run pyKVFinder over a benchmark manifest and normalize residue scores."""

from __future__ import annotations

import argparse
import csv
import json
import shlex
import subprocess
import sys
import time
from pathlib import Path

import residue_score_utils as rsu


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO_ROOT / "data/manifests/p2rank_coach420_mlig_smoke20.csv"
DEFAULT_OUT_DIR = REPO_ROOT / "results/pykvfinder/p2rank_coach420_mlig_smoke20"


PYKVFINDER_CODE = r"""
import argparse
import json
from pathlib import Path

import pyKVFinder


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--json-out", required=True)
    parser.add_argument("--results-out", required=True)
    parser.add_argument("--step", type=float, required=True)
    parser.add_argument("--probe-in", type=float, required=True)
    parser.add_argument("--probe-out", type=float, required=True)
    parser.add_argument("--removal-distance", type=float, required=True)
    parser.add_argument("--volume-cutoff", type=float, required=True)
    parser.add_argument("--ligand-cutoff", type=float, required=True)
    parser.add_argument("--surface", required=True)
    parser.add_argument("--include-depth", action="store_true")
    parser.add_argument("--include-hydropathy", action="store_true")
    parser.add_argument("--ignore-backbone", action="store_true")
    parser.add_argument("--nthreads", type=int)
    args = parser.parse_args()

    result = pyKVFinder.run_workflow(
        input=args.input,
        step=args.step,
        probe_in=args.probe_in,
        probe_out=args.probe_out,
        removal_distance=args.removal_distance,
        volume_cutoff=args.volume_cutoff,
        ligand_cutoff=args.ligand_cutoff,
        include_depth=args.include_depth,
        include_hydropathy=args.include_hydropathy,
        ignore_backbone=args.ignore_backbone,
        surface=args.surface,
        nthreads=args.nthreads,
        verbose=False,
    )
    if result is None:
        payload = {
            "volume": {},
            "area": {},
            "max_depth": {},
            "avg_depth": {},
            "avg_hydropathy": {},
            "residues": {},
            "step": args.step,
            "surface": args.surface,
        }
        Path(args.json_out).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        pyKVFinder.write_results(
            args.results_out,
            input=args.input,
            ligand=None,
            output=None,
            volume={},
            area={},
            max_depth={},
            avg_depth={},
            avg_hydropathy={},
            residues={},
            frequencies=None,
            step=args.step,
        )
        return
    payload = {
        "volume": result.volume or {},
        "area": result.area or {},
        "max_depth": result.max_depth or {},
        "avg_depth": result.avg_depth or {},
        "avg_hydropathy": result.avg_hydropathy or {},
        "residues": result.residues or {},
        "step": args.step,
        "surface": args.surface,
    }
    Path(args.json_out).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    pyKVFinder.write_results(
        args.results_out,
        input=args.input,
        ligand=None,
        output=None,
        volume=result.volume,
        area=result.area,
        max_depth=result.max_depth,
        avg_depth=result.avg_depth,
        avg_hydropathy=result.avg_hydropathy,
        residues=result.residues,
        frequencies=result.frequencies,
        step=args.step,
    )


if __name__ == "__main__":
    main()
"""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Batch-run pyKVFinder cavity detection and write benchmark residue score TSV files.",
    )
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST), help="CSV manifest with target_id and structure_path.")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR), help="Directory for pyKVFinder outputs.")
    parser.add_argument("--python-bin", default=sys.executable, help="Python executable with pyKVFinder installed.")
    parser.add_argument(
        "--score-mode",
        choices=("volume", "depth", "rank-decay", "binary"),
        default="volume",
        help="Cavity score projected to residues. volume/depth are normalized per target.",
    )
    parser.add_argument("--step", type=float, default=0.6, help="pyKVFinder grid step.")
    parser.add_argument("--probe-in", type=float, default=1.4, help="pyKVFinder inner probe radius.")
    parser.add_argument("--probe-out", type=float, default=4.0, help="pyKVFinder outer probe radius.")
    parser.add_argument("--removal-distance", type=float, default=2.4, help="pyKVFinder removal distance.")
    parser.add_argument("--volume-cutoff", type=float, default=5.0, help="Minimum cavity volume retained by pyKVFinder.")
    parser.add_argument("--ligand-cutoff", type=float, default=5.0, help="pyKVFinder ligand cutoff; ligand is not supplied by this adapter.")
    parser.add_argument("--surface", default="SES", choices=("SES", "SAS"), help="pyKVFinder surface mode.")
    parser.add_argument("--include-depth", action="store_true", help="Ask pyKVFinder to calculate cavity depths.")
    parser.add_argument("--include-hydropathy", action="store_true", help="Ask pyKVFinder to calculate hydropathy.")
    parser.add_argument("--ignore-backbone", action="store_true", help="Ignore backbone atoms for pyKVFinder descriptors.")
    parser.add_argument("--nthreads", type=int, help="Thread count passed to pyKVFinder.")
    parser.add_argument("--threshold", type=float, default=0.5, help="Probability threshold recorded in is_binding.")
    parser.add_argument("--limit", type=int, help="Run at most this many manifest rows after filtering.")
    parser.add_argument("--target-id", action="append", default=[], help="Only run this target_id. Can be repeated.")
    parser.add_argument("--timeout", type=float, default=0.0, help="Per-target pyKVFinder timeout in seconds. 0 disables timeout.")
    parser.add_argument("--resume", action="store_true", help="Skip targets whose raw JSON and score TSV already exist.")
    parser.add_argument("--normalize-only", action="store_true", help="Only parse existing pyKVFinder raw JSON outputs.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running pyKVFinder.")
    parser.add_argument("--fail-fast", action="store_true", help="Stop after the first failed target.")
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
    if not rsu.executable_exists(args.python_bin):
        raise SystemExit(f"Python executable not found: {args.python_bin}")
    completed = subprocess.run(
        [args.python_bin, "-c", "import pyKVFinder"],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode:
        raise SystemExit(f"pyKVFinder import failed for {args.python_bin}: {completed.stderr.strip()}")
    if args.score_mode == "depth" and not args.include_depth:
        args.include_depth = True


def run_one(args: argparse.Namespace, row: dict[str, str], out_dir: Path) -> dict[str, str]:
    target_id = rsu.safe_name(row["target_id"])
    input_pdb = rsu.resolve_path(row["structure_path"], REPO_ROOT)
    raw_target_dir = out_dir / "raw" / target_id
    logs_dir = out_dir / "logs"
    scores_tsv = out_dir / "scores" / f"{target_id}.scores.tsv"
    work_pdb = raw_target_dir / f"{target_id}_protein.pdb"
    raw_json = raw_target_dir / "pykvfinder_results.json"
    raw_results = raw_target_dir / "pykvfinder_results.toml"
    stdout_path = logs_dir / f"{target_id}.stdout.txt"
    stderr_path = logs_dir / f"{target_id}.stderr.txt"
    command = build_command(args, work_pdb, raw_json, raw_results)
    command_text = shlex.join(command)
    returncode = 0
    runtime = 0.0
    n_residues = 0

    if args.resume and raw_json.exists() and scores_tsv.exists():
        n_residues = rsu.count_score_rows(scores_tsv)
        return summary(row, input_pdb, raw_target_dir, raw_json, scores_tsv, n_residues, stdout_path, stderr_path, command_text, 0, 0.0, "skipped")

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

        if returncode or not raw_json.exists():
            status = "failed"
        else:
            try:
                n_residues = normalize_residue_scores(
                    input_pdb,
                    raw_json,
                    scores_tsv,
                    row.get("chain_id", "").strip(),
                    args.threshold,
                    args.score_mode,
                )
                status = "ok"
            except ValueError as exc:
                print(f"{target_id}: {exc}", file=sys.stderr)
                status = "failed"

    return summary(
        row,
        input_pdb,
        raw_target_dir,
        raw_json,
        scores_tsv,
        n_residues,
        stdout_path,
        stderr_path,
        command_text,
        returncode,
        runtime,
        status,
    )


def build_command(args: argparse.Namespace, work_pdb: Path, raw_json: Path, raw_results: Path) -> list[str]:
    command = [
        args.python_bin,
        "-c",
        PYKVFINDER_CODE,
        "--input",
        str(work_pdb),
        "--json-out",
        str(raw_json),
        "--results-out",
        str(raw_results),
        "--step",
        str(args.step),
        "--probe-in",
        str(args.probe_in),
        "--probe-out",
        str(args.probe_out),
        "--removal-distance",
        str(args.removal_distance),
        "--volume-cutoff",
        str(args.volume_cutoff),
        "--ligand-cutoff",
        str(args.ligand_cutoff),
        "--surface",
        args.surface,
    ]
    if args.include_depth or args.score_mode == "depth":
        command.append("--include-depth")
    if args.include_hydropathy:
        command.append("--include-hydropathy")
    if args.ignore_backbone:
        command.append("--ignore-backbone")
    if args.nthreads is not None:
        command.extend(["--nthreads", str(args.nthreads)])
    return command


def normalize_residue_scores(
    input_pdb: Path,
    raw_json: Path,
    scores_tsv: Path,
    chain_filter: str,
    threshold: float,
    score_mode: str,
) -> int:
    residue_rows = rsu.load_protein_residues(input_pdb, chain_filter)
    if chain_filter and not residue_rows:
        raise ValueError(f"No protein residues found for chain {chain_filter} in {input_pdb}")

    payload = json.loads(raw_json.read_text(encoding="utf-8"))
    residues_by_cavity = payload.get("residues", {})
    volume_by_cavity = {key: float(value) for key, value in payload.get("volume", {}).items()}
    depth_by_cavity = {key: float(value) for key, value in payload.get("max_depth", {}).items()}
    max_volume = max(volume_by_cavity.values(), default=0.0)
    max_depth = max(depth_by_cavity.values(), default=0.0)

    cavity_ids = sorted(
        residues_by_cavity,
        key=lambda cavity_id: (
            volume_by_cavity.get(cavity_id, 0.0),
            depth_by_cavity.get(cavity_id, 0.0),
            cavity_id,
        ),
        reverse=True,
    )

    residue_scores: dict[tuple[str, str], dict[str, str | float]] = {}
    for rank, cavity_id in enumerate(cavity_ids, start=1):
        probability = cavity_probability(cavity_id, rank, score_mode, volume_by_cavity, depth_by_cavity, max_volume, max_depth)
        for residue_number, chain_id, _resname in residues_by_cavity.get(cavity_id, []):
            key = rsu.residue_key(chain_id, residue_number)
            if chain_filter and key[0] != chain_filter:
                continue
            current = residue_scores.get(key)
            if current is None or probability > float(current["probability"]):
                residue_scores[key] = {"probability": probability, "pocket": cavity_id}

    return rsu.write_residue_scores(scores_tsv, residue_rows, residue_scores, "pykvfinder_score", threshold)


def cavity_probability(
    cavity_id: str,
    rank: int,
    score_mode: str,
    volume_by_cavity: dict[str, float],
    depth_by_cavity: dict[str, float],
    max_volume: float,
    max_depth: float,
) -> float:
    if score_mode == "binary":
        return 1.0
    if score_mode == "rank-decay":
        return 1.0 / rank
    if score_mode == "depth":
        return depth_by_cavity.get(cavity_id, 0.0) / max_depth if max_depth > 0 else 0.0
    return volume_by_cavity.get(cavity_id, 0.0) / max_volume if max_volume > 0 else 0.0


def summary(
    row: dict[str, str],
    input_pdb: Path,
    raw_target_dir: Path,
    raw_json: Path,
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
        "raw_target_dir": str(raw_target_dir),
        "raw_json": str(raw_json),
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
        "raw_target_dir",
        "raw_json",
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
