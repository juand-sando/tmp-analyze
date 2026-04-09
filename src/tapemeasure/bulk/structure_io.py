from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import is_aa

from ..shared.run_logging import get_logger
from .points import normalize_points_dataframe


def _warn(message: str) -> None:
    print(f"WARNING: {message}", file=sys.stderr)
    get_logger().warning(f"warning: {message}")


def _load_structure(input_path: Path):
    suffix = input_path.suffix.lower()
    if suffix == ".pdb":
        parser = PDBParser(QUIET=True)
    elif suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    else:
        raise SystemExit(f"ERROR: unsupported structure file type for {input_path}.")

    try:
        return parser.get_structure(input_path.stem, str(input_path))
    except Exception as exc:
        raise SystemExit(f"ERROR: failed to parse structure file {input_path}: {exc}") from exc


def _select_single_model(structure, input_path: Path):
    models = list(structure.get_models())
    if not models:
        raise SystemExit(f"ERROR: no models found in structure file {input_path}.")
    if len(models) != 1:
        raise SystemExit(
            f"ERROR: structure file {input_path} contains {len(models)} models; "
            "please provide a single-model structure."
        )
    return models[0]


def _protein_chains(model) -> list:
    chains = []
    for chain in model.get_chains():
        if any(is_aa(residue, standard=False) for residue in chain.get_residues()):
            chains.append(chain)
    return chains


def _require_single_chain(model, input_path: Path):
    chains = _protein_chains(model)
    chain_ids = [chain.id.strip() or "<blank>" for chain in chains]
    if not chains:
        raise SystemExit(f"ERROR: no protein chains found in structure file {input_path}.")
    if len(chains) != 1:
        detected = ", ".join(chain_ids)
        _warn(f"structure file {input_path} contains multiple chains: {detected}")
        raise SystemExit(
            f"ERROR: structure file {input_path} must contain exactly one protein chain."
        )
    return chains[0]


def extract_points_from_structure_file(input_path: Path) -> pd.DataFrame:
    structure = _load_structure(input_path)
    model = _select_single_model(structure, input_path)
    chain = _require_single_chain(model, input_path)

    rows: list[dict[str, object]] = []
    missing_ca_count = 0

    for residue in chain.get_residues():
        if not is_aa(residue, standard=False):
            continue
        if "CA" not in residue:
            missing_ca_count += 1
            continue

        atom = residue["CA"]
        coordinates = atom.get_coord()
        insertion_code = residue.id[2].strip() if isinstance(residue.id[2], str) else str(residue.id[2]).strip()
        if insertion_code:
            raise SystemExit(
                f"ERROR: structure file {input_path} contains insertion-coded residues, which are not "
                "supported for the Pos-based smoothing workflow."
            )
        rows.append(
            {
                "Pos": int(residue.id[1]),
                "Res": residue.get_resname().strip(),
                "X": float(coordinates[0]),
                "Y": float(coordinates[1]),
                "Z": float(coordinates[2]),
            }
        )

    if missing_ca_count:
        _warn(
            f"skipped {missing_ca_count} amino-acid residue(s) without a C-alpha atom in {input_path}"
        )

    if not rows:
        raise SystemExit(f"ERROR: no C-alpha coordinates found in structure file {input_path}.")

    return normalize_points_dataframe(
        pd.DataFrame(rows, columns=["Pos", "Res", "X", "Y", "Z"]),
        source_label=str(input_path),
    )
