from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist


def compute_straight_chain_baseline(radius: float, spacing: float = 3.8) -> int:
    if radius <= 0:
        return 0
    expected_neighbors_each_side = int(np.floor(radius / spacing))
    return 2 * expected_neighbors_each_side


def compute_bulkiness_measures(
    points_dataframe: pd.DataFrame,
    radii: list[float],
    k_neighbors: int,
    subtract_baseline: bool,
    baseline_spacing: float,
) -> tuple[pd.DataFrame, dict[str, object]]:
    coordinates = points_dataframe[["X", "Y", "Z"]].to_numpy(dtype=float)
    point_count = len(coordinates)
    if point_count < 2:
        raise SystemExit("ERROR: bulkiness analysis requires at least 2 coordinates.")

    effective_k = min(k_neighbors, point_count - 1)
    kd_tree = cKDTree(coordinates)

    output_dataframe = points_dataframe.copy()
    neighbor_counts_summary: dict[str, dict[str, float | int]] = {}

    for radius in radii:
        neighbors_within_radius = kd_tree.query_ball_point(coordinates, r=radius)
        neighbor_count_within_radius = np.array(
            [len(neighbor_list) - 1 for neighbor_list in neighbors_within_radius],
            dtype=int,
        )

        count_column_name = f"neighbor_count_r{radius:g}"
        output_dataframe[count_column_name] = neighbor_count_within_radius
        neighbor_counts_summary[count_column_name] = {
            "mean": float(np.mean(neighbor_count_within_radius)),
            "median": float(np.median(neighbor_count_within_radius)),
            "min": int(np.min(neighbor_count_within_radius)),
            "max": int(np.max(neighbor_count_within_radius)),
        }

        if subtract_baseline:
            expected_baseline = compute_straight_chain_baseline(radius=radius, spacing=baseline_spacing)
            extra_crowding = neighbor_count_within_radius - expected_baseline
            extra_column_name = f"extra_crowding_r{radius:g}"
            output_dataframe[extra_column_name] = extra_crowding
            neighbor_counts_summary[extra_column_name] = {
                "baseline_expected": int(expected_baseline),
                "mean": float(np.mean(extra_crowding)),
                "median": float(np.median(extra_crowding)),
                "min": int(np.min(extra_crowding)),
                "max": int(np.max(extra_crowding)),
            }

    knn_distances, _ = kd_tree.query(coordinates, k=effective_k + 1)
    first_k_distances = knn_distances[:, 1:effective_k + 1]
    k_mean_neighbor_distance = np.mean(first_k_distances, axis=1)

    output_dataframe["k_mean_neighbor_distance"] = k_mean_neighbor_distance
    output_dataframe["knn_bulkiness_score"] = 1.0 / (k_mean_neighbor_distance ** 3 + 1e-12)

    summary: dict[str, object] = {
        "num_points": int(point_count),
        "radii": radii,
        "k_neighbors_requested": int(k_neighbors),
        "k_neighbors_used": int(effective_k),
        "neighbor_counts": neighbor_counts_summary,
        "k_mean_neighbor_distance": {
            "mean": float(np.mean(k_mean_neighbor_distance)),
            "median": float(np.median(k_mean_neighbor_distance)),
            "min": float(np.min(k_mean_neighbor_distance)),
            "max": float(np.max(k_mean_neighbor_distance)),
        },
    }
    return output_dataframe, summary


def compute_all_pairwise_distances(points_dataframe: pd.DataFrame) -> np.ndarray:
    coordinates = points_dataframe[["X", "Y", "Z"]].to_numpy(dtype=float)
    if len(coordinates) < 2:
        raise SystemExit("ERROR: pairwise distance analysis requires at least 2 coordinates.")
    return pdist(coordinates, metric="euclidean")


def summarize_pairwise_distances(all_pairwise_distances: np.ndarray, point_count: int) -> dict[str, float | int]:
    return {
        "num_points": int(point_count),
        "num_pairs": int(len(all_pairwise_distances)),
        "mean_distance": float(np.mean(all_pairwise_distances)),
        "median_distance": float(np.median(all_pairwise_distances)),
        "min_distance": float(np.min(all_pairwise_distances)),
        "max_distance": float(np.max(all_pairwise_distances)),
    }
