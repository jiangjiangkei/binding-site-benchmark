#!/usr/bin/env python3
"""Run DeepPocket over a benchmark manifest and normalize residue scores."""

from __future__ import annotations

import argparse
import ast
import csv
import os
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
DEFAULT_OUT_DIR = REPO_ROOT / "results/deeppocket/p2rank_coach420_mlig_smoke20"
POCKET_FILE_RE = re.compile(r"pocket(\d+)_atm\.pdb$")
SEGMENTED_POCKET_RE = re.compile(r"_pocket(\d+)\.pdb$")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Batch-run DeepPocket prediction and write benchmark residue score TSV files.",
    )
    parser.add_argument("--manifest", default=str(DEFAULT_MANIFEST), help="CSV manifest with target_id and structure_path.")
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR), help="Directory for DeepPocket outputs.")
    parser.add_argument("--deeppocket-dir", help="Directory containing DeepPocket predict.py.")
    parser.add_argument("--predict-script", help="Explicit path to DeepPocket predict.py.")
    parser.add_argument("--python-bin", default=sys.executable, help="Python executable in the DeepPocket environment.")
    parser.add_argument("--class-checkpoint", help="DeepPocket classification checkpoint.")
    parser.add_argument("--seg-checkpoint", help="DeepPocket segmentation checkpoint.")
    parser.add_argument("--fpocket-bin", help="fpocket executable made available as fpocket on PATH.")
    parser.add_argument("--docker-image", help="Run DeepPocket through this Docker image instead of --python-bin.")
    parser.add_argument("--docker-platform", default="linux/amd64", help="Docker platform used with --docker-image.")
    parser.add_argument("--docker-sudo", action="store_true", help="Prefix Docker commands with sudo -n.")
    parser.add_argument("--docker-deeppocket-dir", default="/opt/DeepPocket", help="DeepPocket directory inside --docker-image.")
    parser.add_argument("--rank", type=int, default=3, help="Number of pockets passed to DeepPocket segmentation.")
    parser.add_argument("--threshold", type=float, default=0.5, help="Residue probability threshold recorded in is_binding.")
    parser.add_argument("--segmentation-threshold", type=float, default=0.5, help="DeepPocket segmentation threshold.")
    parser.add_argument("--mask-dist", type=float, default=3.5, help="DeepPocket mask-to-residue distance.")
    parser.add_argument(
        "--projection-source",
        choices=("ranked-candidates", "segmented-pockets"),
        default="ranked-candidates",
        help="Raw DeepPocket output projected to residue scores.",
    )
    parser.add_argument("--limit", type=int, help="Run at most this many manifest rows after filtering.")
    parser.add_argument("--target-id", action="append", default=[], help="Only run this target_id. Can be repeated.")
    parser.add_argument("--timeout", type=float, default=0.0, help="Per-target DeepPocket timeout in seconds. 0 disables timeout.")
    parser.add_argument("--resume", action="store_true", help="Skip targets whose raw output and score TSV already exist.")
    parser.add_argument("--normalize-only", action="store_true", help="Only parse existing DeepPocket raw outputs.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running DeepPocket.")
    parser.add_argument("--fail-fast", action="store_true", help="Stop after the first failed target.")
    parser.add_argument("--extra-arg", action="append", default=[], help="Extra argument passed through to DeepPocket predict.py.")
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
        missing = []
        for label, value in [
            ("classification checkpoint", args.class_checkpoint),
            ("segmentation checkpoint", args.seg_checkpoint),
        ]:
            if not value or not rsu.resolve_path(value, Path.cwd()).exists():
                missing.append(label)
        if missing:
            raise SystemExit(f"Missing required DeepPocket runtime paths: {', '.join(missing)}")
        return
    predict_script = deep_pocket_predict_script(args)
    missing = []
    for label, value in [
        ("DeepPocket predict.py", predict_script),
        ("classification checkpoint", args.class_checkpoint),
        ("segmentation checkpoint", args.seg_checkpoint),
    ]:
        if not value or not rsu.resolve_path(value, Path.cwd()).exists():
            missing.append(label)
    if missing:
        raise SystemExit(f"Missing required DeepPocket runtime paths: {', '.join(missing)}")
    if not rsu.executable_exists(args.python_bin):
        raise SystemExit(f"Python executable not found: {args.python_bin}")
    if args.fpocket_bin and not rsu.executable_exists(args.fpocket_bin):
        raise SystemExit(f"fpocket executable not found: {args.fpocket_bin}")


def run_one(args: argparse.Namespace, row: dict[str, str], out_dir: Path) -> dict[str, str]:
    target_id = rsu.safe_name(row["target_id"])
    input_pdb = rsu.resolve_path(row["structure_path"], REPO_ROOT)
    raw_target_dir = out_dir / "raw" / target_id
    logs_dir = out_dir / "logs"
    scores_tsv = out_dir / "scores" / f"{target_id}.scores.tsv"
    work_pdb = raw_target_dir / f"{target_id}.pdb"
    pockets_dir = raw_target_dir / f"{target_id}_nowat_out" / "pockets"
    stdout_path = logs_dir / f"{target_id}.stdout.txt"
    stderr_path = logs_dir / f"{target_id}.stderr.txt"
    command = build_command(args, work_pdb)
    command_text = shlex.join(command)
    returncode = 0
    runtime = 0.0
    n_residues = 0

    if args.resume and raw_output_exists(raw_target_dir, pockets_dir) and scores_tsv.exists():
        n_residues = rsu.count_score_rows(scores_tsv)
        return summary(row, input_pdb, raw_target_dir, pockets_dir, scores_tsv, n_residues, stdout_path, stderr_path, command_text, 0, 0.0, "skipped")

    raw_target_dir.mkdir(parents=True, exist_ok=True)
    if not args.dry_run:
        shutil.copy2(input_pdb, work_pdb)

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
                    env=runtime_env(args),
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
                    raw_target_dir,
                    pockets_dir,
                    scores_tsv,
                    row.get("chain_id", "").strip(),
                    args.threshold,
                    args.projection_source,
                )
                status = "ok"
            except ValueError as exc:
                print(f"{target_id}: {exc}", file=sys.stderr)
                status = "failed"

    return summary(
        row,
        input_pdb,
        raw_target_dir,
        pockets_dir,
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
                f"{Path.home()}:{Path.home()}",
                "-v",
                "/tmp:/tmp",
                "-w",
                str(work_pdb.parent),
                args.docker_image,
                str(Path(args.docker_deeppocket_dir) / "predict.py"),
                "-p",
                str(work_pdb),
                "-c",
                str(rsu.resolve_path(args.class_checkpoint, Path.cwd())) if args.class_checkpoint else "",
                "-s",
                str(rsu.resolve_path(args.seg_checkpoint, Path.cwd())) if args.seg_checkpoint else "",
                "-r",
                str(args.rank),
                "-t",
                str(args.segmentation_threshold),
                "--mask_dist",
                str(args.mask_dist),
                *args.extra_arg,
            ]
        )
        return command

    return [
        args.python_bin,
        str(deep_pocket_predict_script(args)),
        "-p",
        str(work_pdb),
        "-c",
        str(rsu.resolve_path(args.class_checkpoint, Path.cwd())) if args.class_checkpoint else "",
        "-s",
        str(rsu.resolve_path(args.seg_checkpoint, Path.cwd())) if args.seg_checkpoint else "",
        "-r",
        str(args.rank),
        "-t",
        str(args.segmentation_threshold),
        "--mask_dist",
        str(args.mask_dist),
        *args.extra_arg,
    ]


def deep_pocket_predict_script(args: argparse.Namespace) -> Path:
    if args.predict_script:
        return rsu.resolve_path(args.predict_script, Path.cwd())
    if not args.deeppocket_dir:
        return Path("predict.py")
    return rsu.resolve_path(args.deeppocket_dir, Path.cwd()) / "predict.py"


def runtime_env(args: argparse.Namespace) -> dict[str, str]:
    env = os.environ.copy()
    if args.fpocket_bin:
        fpocket_path = rsu.resolve_path(args.fpocket_bin, Path.cwd())
        env["PATH"] = str(fpocket_path.parent) + os.pathsep + env.get("PATH", "")
    return env


def normalize_residue_scores(
    input_pdb: Path,
    raw_target_dir: Path,
    pockets_dir: Path,
    scores_tsv: Path,
    chain_filter: str,
    threshold: float,
    projection_source: str,
) -> int:
    residue_rows = rsu.load_protein_residues(input_pdb, chain_filter)
    if chain_filter and not residue_rows:
        raise ValueError(f"No protein residues found for chain {chain_filter} in {input_pdb}")

    if projection_source == "segmented-pockets":
        residue_scores = segmented_pocket_scores(raw_target_dir, chain_filter)
    else:
        residue_scores = ranked_candidate_scores(pockets_dir, chain_filter)

    return rsu.write_residue_scores(scores_tsv, residue_rows, residue_scores, "deeppocket_score", threshold)


def ranked_candidate_scores(pockets_dir: Path, chain_filter: str) -> dict[tuple[str, str], dict[str, str | float]]:
    ranked_types = pockets_dir / "bary_centers_ranked.types"
    confidence_path = pockets_dir / "bary_centers_confidence.txt"
    if not ranked_types.exists():
        raise ValueError(f"Missing DeepPocket ranked types file: {ranked_types}")

    confidences = parse_confidences(confidence_path)
    residue_scores: dict[tuple[str, str], dict[str, str | float]] = {}
    with ranked_types.open(encoding="utf-8", errors="replace") as handle:
        for rank_index, line in enumerate(handle, start=1):
            parts = line.split()
            if not parts:
                continue
            pocket_index = parts[0]
            probability = confidences[rank_index - 1] if rank_index - 1 < len(confidences) else 1.0 / rank_index
            pocket_path = pockets_dir / f"pocket{pocket_index}_atm.pdb"
            if not pocket_path.exists():
                continue
            for residue in rsu.residues_from_pdb(pocket_path, chain_filter):
                current = residue_scores.get(residue)
                if current is None or probability > float(current["probability"]):
                    residue_scores[residue] = {"probability": probability, "pocket": pocket_index}
    return residue_scores


def segmented_pocket_scores(raw_target_dir: Path, chain_filter: str) -> dict[tuple[str, str], dict[str, str | float]]:
    pockets = sorted(raw_target_dir.glob("*_pocket*.pdb"), key=segmented_sort_key)
    residue_scores: dict[tuple[str, str], dict[str, str | float]] = {}
    for rank_index, pocket_path in enumerate(pockets, start=1):
        probability = 1.0 / rank_index
        pocket = str(rank_index)
        for residue in rsu.residues_from_pdb(pocket_path, chain_filter):
            current = residue_scores.get(residue)
            if current is None or probability > float(current["probability"]):
                residue_scores[residue] = {"probability": probability, "pocket": pocket}
    return residue_scores


def parse_confidences(path: Path) -> list[float]:
    if not path.exists():
        return []
    text = path.read_text(encoding="utf-8", errors="replace").strip()
    if not text:
        return []
    try:
        value = ast.literal_eval(text)
        return [float(item) for item in value]
    except (SyntaxError, ValueError, TypeError):
        return [float(match) for match in re.findall(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?", text)]


def raw_output_exists(raw_target_dir: Path, pockets_dir: Path) -> bool:
    return (pockets_dir / "bary_centers_ranked.types").exists() or any(raw_target_dir.glob("*_pocket*.pdb"))


def segmented_sort_key(path: Path) -> int:
    match = SEGMENTED_POCKET_RE.search(path.name)
    return int(match.group(1)) if match else 10**9


def summary(
    row: dict[str, str],
    input_pdb: Path,
    raw_target_dir: Path,
    pockets_dir: Path,
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
        "pockets_dir": str(pockets_dir),
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
        "pockets_dir",
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
