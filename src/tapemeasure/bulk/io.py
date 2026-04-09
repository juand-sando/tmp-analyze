from __future__ import annotations

from pathlib import Path

import pandas as pd

from .structure_io import extract_points_from_structure_file

REQUIRED_COLUMNS = {"Pos", "Res", "X", "Y", "Z"}


def load_points(csv_path: Path) -> pd.DataFrame:
    points_dataframe = pd.read_csv(csv_path)
    if not REQUIRED_COLUMNS.issubset(points_dataframe.columns):
        raise SystemExit(
            f"ERROR: {csv_path} must contain columns {sorted(REQUIRED_COLUMNS)}, "
            f"but has {sorted(points_dataframe.columns)}."
        )
    return points_dataframe


def load_points_from_input(input_path: Path) -> tuple[pd.DataFrame, str]:
    suffix = input_path.suffix.lower()
    if suffix in {".csv"}:
        return load_points(input_path), "csv"
    if suffix in {".pdb", ".cif", ".mmcif"}:
        return extract_points_from_structure_file(input_path), "structure"
    raise SystemExit(
        f"ERROR: unsupported input file type for {input_path}. "
        "Expected .csv, .pdb, .cif, or .mmcif."
    )
