"""Shared helpers for projecting pocket-style outputs to residue score TSVs."""

from __future__ import annotations

import csv
import math
import os
import re
import shutil
from pathlib import Path


PDB_ATOM_RECORDS = {"ATOM", "HETATM"}


def load_manifest(path: Path, selected_target_ids: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"target_id", "structure_path"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Manifest missing required columns: {', '.join(sorted(missing))}")
        return [row for row in reader if not selected_target_ids or row["target_id"] in selected_target_ids]


def load_protein_residues(path: Path, chain_filter: str) -> list[dict[str, str]]:
    seen = set()
    rows = []
    for atom in load_protein_atoms(path, chain_filter):
        key = residue_key(atom["chain_id"], atom["residue_number"])
        if key in seen:
            continue
        seen.add(key)
        rows.append({"chain_id": key[0], "residue_number": key[1]})
    return rows


def load_protein_atoms(path: Path, chain_filter: str) -> list[dict[str, object]]:
    atoms = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line[:6].strip() != "ATOM":
                continue
            atom = parse_pdb_atom(line)
            if atom is None:
                continue
            if chain_filter and atom["chain_id"] != chain_filter:
                continue
            atoms.append(atom)
    return atoms


def parse_pdb_atom(line: str) -> dict[str, object] | None:
    if len(line) < 54:
        return None
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
    except ValueError:
        return None

    residue_number = line[22:26].strip()
    insertion_code = line[26].strip()
    if insertion_code:
        residue_number = f"{residue_number}{insertion_code}"
    atom_name = line[12:16].strip()
    element = line[76:78].strip() or "".join(char for char in atom_name if char.isalpha())[:1]
    return {
        "record": line[:6].strip(),
        "atom_name": atom_name,
        "element": element.upper(),
        "chain_id": line[21].strip(),
        "resname": line[17:20].strip(),
        "residue_number": residue_number,
        "xyz": (x, y, z),
    }


def residues_from_pdb(path: Path, chain_filter: str) -> set[tuple[str, str]]:
    residues = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line[:6].strip() not in PDB_ATOM_RECORDS:
                continue
            atom = parse_pdb_atom(line)
            if atom is None:
                continue
            if chain_filter and atom["chain_id"] != chain_filter:
                continue
            residues.add(residue_key(atom["chain_id"], atom["residue_number"]))
    return residues


def residue_scores_from_pdb_atoms(path: Path, chain_filter: str) -> dict[tuple[str, str], float]:
    scores: dict[tuple[str, str], float] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line[:6].strip() not in PDB_ATOM_RECORDS:
                continue
            atom = parse_pdb_atom(line)
            if atom is None:
                continue
            if chain_filter and atom["chain_id"] != chain_filter:
                continue
            try:
                score = float(line[54:60])
            except ValueError:
                score = 1.0
            key = residue_key(atom["chain_id"], atom["residue_number"])
            scores[key] = max(scores.get(key, 0.0), score)
    return scores


def mol2_atom_coordinates(path: Path) -> list[tuple[float, float, float]]:
    coordinates = []
    in_atoms = False
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped.upper() == "@<TRIPOS>ATOM":
                in_atoms = True
                continue
            if stripped.startswith("@<TRIPOS>") and in_atoms:
                break
            if not in_atoms or not stripped:
                continue
            parts = stripped.split()
            if len(parts) < 5:
                continue
            try:
                coordinates.append((float(parts[2]), float(parts[3]), float(parts[4])))
            except ValueError:
                continue
    return coordinates


def residues_from_coordinates(
    coordinates: list[tuple[float, float, float]],
    protein_atoms: list[dict[str, object]],
    max_distance: float,
) -> set[tuple[str, str]]:
    residues = set()
    if not coordinates or not protein_atoms:
        return residues

    exact_index: dict[tuple[int, int, int], list[dict[str, object]]] = {}
    for atom in protein_atoms:
        exact_index.setdefault(rounded_xyz(atom["xyz"], 3), []).append(atom)

    max_distance2 = max_distance * max_distance
    for xyz in coordinates:
        candidates = exact_index.get(rounded_xyz(xyz, 3), [])
        if not candidates:
            candidates = protein_atoms
        best_atom = None
        best_distance2 = math.inf
        x, y, z = xyz
        for atom in candidates:
            ax, ay, az = atom["xyz"]
            distance2 = (x - ax) ** 2 + (y - ay) ** 2 + (z - az) ** 2
            if distance2 < best_distance2:
                best_distance2 = distance2
                best_atom = atom
        if best_atom is not None and best_distance2 <= max_distance2:
            residues.add(residue_key(best_atom["chain_id"], best_atom["residue_number"]))
    return residues


def write_residue_scores(
    scores_tsv: Path,
    residue_rows: list[dict[str, str]],
    residue_scores: dict[tuple[str, str], dict[str, str | float]],
    method_score_field: str,
    threshold: float,
) -> int:
    scores_tsv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["residue_id", "chain_id", "residue_number", "probability", "is_binding", method_score_field, "pocket"]
    rows = []
    for residue in residue_rows:
        key = residue_key(residue["chain_id"], residue["residue_number"])
        score_row = residue_scores.get(key, {})
        probability = float(score_row.get("probability", 0.0))
        pocket = str(score_row.get("pocket", "0"))
        rows.append(
            {
                "residue_id": f"{residue['chain_id']}_{residue['residue_number']}",
                "chain_id": residue["chain_id"],
                "residue_number": residue["residue_number"],
                "probability": f"{probability:.6f}",
                "is_binding": "1" if probability >= threshold else "0",
                method_score_field: f"{probability:.6f}",
                "pocket": pocket,
            }
        )

    with scores_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def write_protein_only_pdb(input_pdb: Path, output_pdb: Path, chain_filter: str) -> None:
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    wrote_any = False
    with input_pdb.open(encoding="utf-8", errors="replace") as source, output_pdb.open("w", encoding="utf-8") as dest:
        for line in source:
            if line[:6].strip() != "ATOM":
                continue
            atom = parse_pdb_atom(line)
            if atom is None:
                continue
            if chain_filter and atom["chain_id"] != chain_filter:
                continue
            dest.write(line.rstrip("\n") + "\n")
            wrote_any = True
        if wrote_any:
            dest.write("END\n")


def count_score_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open(newline="", encoding="utf-8") as handle:
        return max(0, sum(1 for _ in handle) - 1)


def rounded_xyz(xyz: object, digits: int) -> tuple[int, int, int]:
    x, y, z = xyz
    scale = 10**digits
    return round(float(x) * scale), round(float(y) * scale), round(float(z) * scale)


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
