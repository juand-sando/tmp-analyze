from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from ..shared.run_logging import log_file_write


def _configure_matplotlib():
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    mpl.rcParams["svg.fonttype"] = "none"
    return plt


def plot_bulkiness(
    output_dataframe: pd.DataFrame,
    output_dir: Path,
    file_prefix: str = "",
) -> list[Path]:
    plt = _configure_matplotlib()
    written_paths: list[Path] = []
    prefix = f"{file_prefix}_" if file_prefix else ""

    bulkiness_score_path = output_dir / f"{prefix}bulkiness_score_profile.png"
    position_values = output_dataframe["Pos"].to_numpy(dtype=int)
    if position_values.size == 0:
        return written_paths

    full_position_range = np.arange(int(position_values.min()), int(position_values.max()) + 1, dtype=int)
    bulkiness_series = (
        output_dataframe.set_index("Pos")["bulkiness_score_ca_per_nm3"]
        .reindex(full_position_range)
    )

    plt.figure()
    plt.plot(full_position_range, bulkiness_series.to_numpy(dtype=float), color="#0b5fa5", linewidth=1.6)
    plt.xlabel("Pos (residue index)")
    plt.ylabel("Bulkiness score (Ca/nm^3)")
    plt.title("Bulkiness score along chain")
    plt.tight_layout()
    log_file_write(bulkiness_score_path)
    plt.savefig(bulkiness_score_path, dpi=200)
    plt.close()
    written_paths.append(bulkiness_score_path)

    return written_paths


def plot_nearest_neighbor_profiles(
    profile_dataframe: pd.DataFrame,
    max_neighbors: int,
    output_png_path: Path,
) -> None:
    plt = _configure_matplotlib()
    point_index = profile_dataframe["Pos"].to_numpy()
    column_count = 3
    row_count = int(math.ceil(max_neighbors / column_count))
    figure, axes = plt.subplots(row_count, column_count, figsize=(15, 3.2 * row_count), sharex=True)
    axes_array = np.atleast_1d(axes).ravel()

    for neighbor_rank in range(1, max_neighbors + 1):
        axis = axes_array[neighbor_rank - 1]
        raw_column_name = f"neighbor_{neighbor_rank}_distance"
        axis.plot(point_index, profile_dataframe[raw_column_name], color="#0b5fa5", linewidth=1.3)
        axis.set_title(f"{neighbor_rank}th neighbor", fontsize=10)
        axis.set_ylabel("Distance")
        axis.grid(alpha=0.2, linewidth=0.5)

    for unused_axis in axes_array[max_neighbors:]:
        unused_axis.axis("off")

    for axis in axes_array[-column_count:]:
        if axis.has_data():
            axis.set_xlabel("Pos (residue index)")

    figure.suptitle("Nearest-neighbor distance profiles", fontsize=14)
    figure.tight_layout(rect=(0, 0, 1, 0.98))
    log_file_write(output_png_path)
    figure.savefig(output_png_path, dpi=200)
    plt.close(figure)
