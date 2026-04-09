from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def sphere_volume_nm3(radius_angstrom: float) -> float:
    radius_nm = radius_angstrom / 10.0
    return (4.0 / 3.0) * np.pi * (radius_nm ** 3)


def compute_bulkiness_measures(
    points_dataframe: pd.DataFrame,
    radius: float,
) -> tuple[pd.DataFrame, dict[str, object]]:
    coordinates = points_dataframe[["X", "Y", "Z"]].to_numpy(dtype=float)
    point_count = len(coordinates)
    if point_count < 2:
        raise SystemExit("ERROR: bulkiness analysis requires at least 2 coordinates.")

    kd_tree = cKDTree(coordinates)
    output_dataframe = points_dataframe.copy()
    neighbors_within_radius = kd_tree.query_ball_point(coordinates, r=radius)
    neighbor_count_within_radius = np.array(
        [len(neighbor_list) - 1 for neighbor_list in neighbors_within_radius],
        dtype=int,
    )

    output_dataframe["neighbor_count"] = neighbor_count_within_radius
    smoothed_values = iterative_chain_smoothing(output_dataframe["Pos"].tolist(), neighbor_count_within_radius, rounds=3)
    output_dataframe["neighbor_count_smoothed"] = output_dataframe["Pos"].map(smoothed_values)
    bulkiness_denominator_nm3 = sphere_volume_nm3(radius)
    output_dataframe["bulkiness_score_ca_per_nm3"] = (
        output_dataframe["neighbor_count_smoothed"] / bulkiness_denominator_nm3
    )

    smoothed_array = output_dataframe["neighbor_count_smoothed"].dropna().to_numpy(dtype=float)
    bulkiness_score_array = output_dataframe["bulkiness_score_ca_per_nm3"].dropna().to_numpy(dtype=float)

    summary: dict[str, object] = {
        "num_points": int(point_count),
        "radius": float(radius),
        "radius_nm": float(radius / 10.0),
        "sphere_volume_nm3": float(bulkiness_denominator_nm3),
        "smoothing_rounds": 3,
        "neighbor_count": {
            "mean": float(np.mean(neighbor_count_within_radius)),
            "median": float(np.median(neighbor_count_within_radius)),
            "min": int(np.min(neighbor_count_within_radius)),
            "max": int(np.max(neighbor_count_within_radius)),
        },
        "smoothed_neighbor_count": None
        if smoothed_array.size == 0
        else {
            "mean": float(np.mean(smoothed_array)),
            "median": float(np.median(smoothed_array)),
            "min": float(np.min(smoothed_array)),
            "max": float(np.max(smoothed_array)),
            "surviving_positions": int(smoothed_array.size),
        },
        "bulkiness_score_ca_per_nm3": None
        if bulkiness_score_array.size == 0
        else {
            "mean": float(np.mean(bulkiness_score_array)),
            "median": float(np.median(bulkiness_score_array)),
            "min": float(np.min(bulkiness_score_array)),
            "max": float(np.max(bulkiness_score_array)),
            "surviving_positions": int(bulkiness_score_array.size),
        },
    }
    return output_dataframe, summary


def iterative_chain_smoothing(
    positions: list[int],
    values: np.ndarray,
    *,
    rounds: int,
) -> dict[int, float]:
    current_values = {int(position): float(value) for position, value in zip(positions, values, strict=True)}
    for _ in range(rounds):
        next_values: dict[int, float] = {}
        for position in sorted(current_values.keys()):
            left_position = position - 1
            right_position = position + 1
            if left_position not in current_values or right_position not in current_values:
                continue
            next_values[position] = (
                0.5 * current_values[position]
                + 0.25 * current_values[left_position]
                + 0.25 * current_values[right_position]
            )
        current_values = next_values
        if not current_values:
            break
    return current_values


def compute_nearest_neighbor_distance_profiles(
    points_dataframe: pd.DataFrame,
    max_neighbors: int,
) -> tuple[pd.DataFrame, dict[str, object]]:
    coordinates = points_dataframe[["X", "Y", "Z"]].to_numpy(dtype=float)
    point_count = len(coordinates)
    if point_count < 2:
        raise SystemExit("ERROR: nearest-neighbor analysis requires at least 2 coordinates.")

    effective_max_neighbors = min(max_neighbors, point_count - 1)
    kd_tree = cKDTree(coordinates)
    neighbor_distances, _ = kd_tree.query(coordinates, k=effective_max_neighbors + 1)
    ranked_neighbor_distances = neighbor_distances[:, 1:effective_max_neighbors + 1]

    output_dataframe = points_dataframe[["Pos", "Res"]].copy()
    summary_by_rank: dict[str, dict[str, float | int]] = {}

    for neighbor_rank in range(1, effective_max_neighbors + 1):
        distance_values = ranked_neighbor_distances[:, neighbor_rank - 1]
        output_dataframe[f"neighbor_{neighbor_rank}_distance"] = distance_values
        summary_by_rank[f"neighbor_{neighbor_rank}"] = {
            "mean_distance": float(np.mean(distance_values)),
            "median_distance": float(np.median(distance_values)),
            "min_distance": float(np.min(distance_values)),
            "max_distance": float(np.max(distance_values)),
        }

    summary: dict[str, object] = {
        "num_points": int(point_count),
        "max_neighbors_requested": int(max_neighbors),
        "max_neighbors_used": int(effective_max_neighbors),
        "neighbor_distance_summary": summary_by_rank,
    }
    return output_dataframe, summary
