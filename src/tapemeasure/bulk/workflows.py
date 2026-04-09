from __future__ import annotations

import pandas as pd

from .analysis import compute_all_pairwise_distances, compute_bulkiness_measures, summarize_pairwise_distances
from .io import load_points
from .models import BulkinessConfig
from .outputs import analysis_dir, write_analysis_metadata, write_dataframe_csv, write_summary_json
from .plotting import plot_bulkiness, plot_pairwise_distance_histogram


def run(config: BulkinessConfig) -> str:
    write_analysis_metadata(config)
    if config.mode == "bulkiness":
        return run_bulkiness(config)
    if config.mode == "pairwise-histogram":
        return run_pairwise_histogram(config)
    raise SystemExit(f"ERROR: unsupported mode {config.mode}")


def run_bulkiness(config: BulkinessConfig) -> str:
    points_dataframe = load_points(config.input_csv)
    output_dataframe, summary = compute_bulkiness_measures(
        points_dataframe=points_dataframe,
        radii=config.radii,
        k_neighbors=config.k_neighbors,
        subtract_baseline=config.subtract_baseline,
        baseline_spacing=config.baseline_spacing,
    )

    output_directory = analysis_dir(config.output_dir)
    output_directory.mkdir(parents=True, exist_ok=True)
    write_dataframe_csv(output_dataframe, output_directory / "bulkiness_table.csv")
    write_summary_json(summary, output_directory / "bulkiness_summary.json")

    plot_count = 0
    if config.make_plots:
        plot_count = len(
            plot_bulkiness(
                output_dataframe=output_dataframe,
                radii=config.radii,
                subtract_baseline=config.subtract_baseline,
                output_dir=output_directory,
            )
        )

    return (
        f"Done. Processed {len(points_dataframe)} coordinate(s); "
        f"wrote bulkiness table and summary; "
        f"wrote {plot_count} plot(s)."
    )


def run_pairwise_histogram(config: BulkinessConfig) -> str:
    points_dataframe = load_points(config.input_csv)
    all_pairwise_distances = compute_all_pairwise_distances(points_dataframe)
    summary = summarize_pairwise_distances(all_pairwise_distances, len(points_dataframe))

    output_directory = analysis_dir(config.output_dir)
    output_directory.mkdir(parents=True, exist_ok=True)
    write_dataframe_csv(
        pd.DataFrame({"pairwise_distance": all_pairwise_distances}),
        output_directory / "pairwise_distances.csv",
    )
    write_summary_json(summary, output_directory / "pairwise_distance_summary.json")
    plot_pairwise_distance_histogram(
        all_pairwise_distances=all_pairwise_distances,
        histogram_bins=config.histogram_bins,
        output_png_path=output_directory / "pairwise_distance_hist.png",
    )

    return (
        f"Done. Processed {len(points_dataframe)} coordinate(s); "
        f"computed {len(all_pairwise_distances)} pairwise distance(s); "
        f"wrote histogram outputs."
    )
