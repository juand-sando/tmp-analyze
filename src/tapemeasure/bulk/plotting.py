from __future__ import annotations

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
    radii: list[float],
    subtract_baseline: bool,
    output_dir: Path,
    file_prefix: str = "",
) -> list[Path]:
    plt = _configure_matplotlib()
    point_index = output_dataframe["Pos"].to_numpy()
    written_paths: list[Path] = []
    prefix = f"{file_prefix}_" if file_prefix else ""

    neighbor_counts_path = output_dir / f"{prefix}neighbor_counts.png"
    plt.figure()
    for radius in radii:
        plt.plot(point_index, output_dataframe[f"neighbor_count_r{radius:g}"], label=f"count r={radius:g}")
    plt.xlabel("Pos (residue index)")
    plt.ylabel("Neighbors within radius")
    plt.title("Local neighbor counts along chain")
    plt.legend()
    plt.tight_layout()
    log_file_write(neighbor_counts_path)
    plt.savefig(neighbor_counts_path, dpi=200)
    plt.close()
    written_paths.append(neighbor_counts_path)

    if subtract_baseline:
        extra_crowding_path = output_dir / f"{prefix}extra_crowding.png"
        plt.figure()
        for radius in radii:
            plt.plot(point_index, output_dataframe[f"extra_crowding_r{radius:g}"], label=f"extra r={radius:g}")
        plt.xlabel("Pos (residue index)")
        plt.ylabel("Extra neighbors above straight-chain baseline")
        plt.title("Extra crowding along chain")
        plt.legend()
        plt.tight_layout()
        log_file_write(extra_crowding_path)
        plt.savefig(extra_crowding_path, dpi=200)
        plt.close()
        written_paths.append(extra_crowding_path)

    knn_bulkiness_path = output_dir / f"{prefix}knn_bulkiness.png"
    plt.figure()
    plt.plot(point_index, output_dataframe["knn_bulkiness_score"])
    plt.xlabel("Pos (residue index)")
    plt.ylabel("Bulkiness ~ 1 / (mean kNN distance)^3")
    plt.title("kNN bulkiness along chain")
    plt.tight_layout()
    log_file_write(knn_bulkiness_path)
    plt.savefig(knn_bulkiness_path, dpi=200)
    plt.close()
    written_paths.append(knn_bulkiness_path)

    k_mean_distance_hist_path = output_dir / f"{prefix}k_mean_distance_hist.png"
    plt.figure()
    plt.hist(output_dataframe["k_mean_neighbor_distance"], bins=30)
    plt.xlabel("Mean distance to first k neighbors")
    plt.ylabel("Count")
    plt.title("Distribution of k-mean neighbor distances")
    plt.tight_layout()
    log_file_write(k_mean_distance_hist_path)
    plt.savefig(k_mean_distance_hist_path, dpi=200)
    plt.close()
    written_paths.append(k_mean_distance_hist_path)

    return written_paths


def plot_pairwise_distance_histogram(
    all_pairwise_distances: np.ndarray,
    histogram_bins: int,
    output_png_path: Path,
) -> None:
    plt = _configure_matplotlib()
    plt.figure()
    plt.hist(all_pairwise_distances, bins=histogram_bins)
    plt.xlabel("Pairwise distance")
    plt.ylabel("Count")
    plt.title("Histogram of all pairwise C-alpha distances")
    plt.tight_layout()
    log_file_write(output_png_path)
    plt.savefig(output_png_path, dpi=200)
    plt.close()
