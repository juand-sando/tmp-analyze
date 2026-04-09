from __future__ import annotations

import pandas as pd

from .analysis import compute_bulkiness_measures, compute_nearest_neighbor_distance_profiles
from .io import load_points_from_input
from .models import BulkinessConfig
from .outputs import analysis_dir, write_analysis_metadata, write_dataframe_csv, write_summary_json
from .plotting import plot_bulkiness, plot_nearest_neighbor_profiles


def run(config: BulkinessConfig) -> str:
    write_analysis_metadata(config)
    if config.mode == "bulkiness":
        return run_bulkiness(config)
    if config.mode == "pairwise-histogram":
        return run_pairwise_histogram(config)
    raise SystemExit(f"ERROR: unsupported mode {config.mode}")


def run_bulkiness(config: BulkinessConfig) -> str:
    output_directory = analysis_dir(config.output_dir)
    output_directory.mkdir(parents=True, exist_ok=True)

    points_dataframe, input_kind = load_points_from_input(config.input_path)
    if input_kind == "structure":
        write_dataframe_csv(points_dataframe, output_directory / "extracted_ca_coordinates.csv")

    output_dataframe, summary = compute_bulkiness_measures(
        points_dataframe=points_dataframe,
        radius=config.radius,
    )

    write_dataframe_csv(output_dataframe, output_directory / "bulkiness_table.csv")
    write_summary_json(summary, output_directory / "bulkiness_summary.json")

    plot_count = 0
    if config.make_plots:
        plot_count = len(
            plot_bulkiness(
                output_dataframe=output_dataframe,
                output_dir=output_directory,
            )
        )

    return (
        f"Done. Processed {len(points_dataframe)} coordinate(s) from {input_kind} input; "
        f"wrote bulkiness table and summary; "
        f"wrote {plot_count} plot(s)."
    )


def run_pairwise_histogram(config: BulkinessConfig) -> str:
    output_directory = analysis_dir(config.output_dir)
    output_directory.mkdir(parents=True, exist_ok=True)

    points_dataframe, input_kind = load_points_from_input(config.input_path)
    if input_kind == "structure":
        write_dataframe_csv(points_dataframe, output_directory / "extracted_ca_coordinates.csv")

    profile_dataframe, summary = compute_nearest_neighbor_distance_profiles(
        points_dataframe=points_dataframe,
        max_neighbors=config.max_neighbors,
    )

    write_dataframe_csv(
        profile_dataframe,
        output_directory / "nearest_neighbor_distance_profiles.csv",
    )
    write_summary_json(summary, output_directory / "nearest_neighbor_distance_summary.json")
    plot_nearest_neighbor_profiles(
        profile_dataframe=profile_dataframe,
        max_neighbors=int(summary["max_neighbors_used"]),
        output_png_path=output_directory / "nearest_neighbor_distance_profiles.png",
    )

    return (
        f"Done. Processed {len(points_dataframe)} coordinate(s) from {input_kind} input; "
        f"computed nearest-neighbor profiles through rank {summary['max_neighbors_used']}; "
        f"wrote profile outputs."
    )
