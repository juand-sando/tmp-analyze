from __future__ import annotations

import pandas as pd


def normalize_points_dataframe(points_dataframe: pd.DataFrame, *, source_label: str) -> pd.DataFrame:
    normalized_dataframe = points_dataframe.copy()
    try:
        normalized_dataframe["Pos"] = pd.to_numeric(normalized_dataframe["Pos"], errors="raise").astype(int)
    except Exception as exc:
        raise SystemExit(f"ERROR: {source_label} contains non-integer Pos values.") from exc

    duplicate_positions = normalized_dataframe["Pos"][normalized_dataframe["Pos"].duplicated()].unique().tolist()
    if duplicate_positions:
        duplicate_text = ", ".join(str(value) for value in sorted(duplicate_positions))
        raise SystemExit(f"ERROR: {source_label} contains duplicate Pos values: {duplicate_text}")

    normalized_dataframe = normalized_dataframe.sort_values("Pos").reset_index(drop=True)
    return normalized_dataframe
