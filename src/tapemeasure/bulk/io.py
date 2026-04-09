from __future__ import annotations

from pathlib import Path

import pandas as pd


REQUIRED_COLUMNS = {"Pos", "Res", "X", "Y", "Z"}


def load_points(csv_path: Path) -> pd.DataFrame:
    points_dataframe = pd.read_csv(csv_path)
    if not REQUIRED_COLUMNS.issubset(points_dataframe.columns):
        raise SystemExit(
            f"ERROR: {csv_path} must contain columns {sorted(REQUIRED_COLUMNS)}, "
            f"but has {sorted(points_dataframe.columns)}."
        )
    return points_dataframe
