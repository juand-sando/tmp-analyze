from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from .analysis import compute_bulkiness_measures, compute_nearest_neighbor_distance_profiles
from .io import load_points
from .plotting import plot_bulkiness, plot_nearest_neighbor_profiles


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="Compute radius-based neighbor counts or nearest-neighbor profiles from a 3D C-alpha CSV."
    )
    parser.add_argument("csv_path", type=str, help="Path to input CSV with columns Pos, Res, X, Y, Z.")
    parser.add_argument(
        "--mode",
        type=str,
        choices=["bulkiness", "pairwise_histogram"],
        default="bulkiness",
        help="Which output to produce.",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=10,
        help="Radius in Angstroms used to count spatial neighbors.",
    )
    parser.add_argument(
        "--output_csv",
        type=str,
        default=None,
        help="Output CSV path for bulkiness mode. Default: <input>_bulkiness.csv",
    )
    parser.add_argument(
        "--make_plots",
        action="store_true",
        help="Save PNG plots for bulkiness mode.",
    )
    parser.add_argument(
        "--output_prefix",
        type=str,
        default=None,
        help="Prefix for bulkiness plot filenames. Default: output CSV stem.",
    )
    parser.add_argument(
        "--max_neighbors",
        type=int,
        default=15,
        help="Maximum nearest-neighbor rank to analyze.",
    )
    parser.add_argument(
        "--pairwise_output_csv",
        type=str,
        default=None,
        help="Output CSV path for pairwise mode. Default: <input>_nearest_neighbor_distance_profiles.csv",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_arguments(argv)
    csv_path = Path(args.csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Input file not found: {csv_path}")

    points_dataframe = load_points(csv_path)

    if args.mode == "bulkiness":
        output_dataframe, summary = compute_bulkiness_measures(
            points_dataframe=points_dataframe,
            radius=args.radius,
        )

        output_csv_path = Path(args.output_csv) if args.output_csv else csv_path.with_name(csv_path.stem + "_bulkiness.csv")
        output_dataframe.to_csv(output_csv_path, index=False)

        print("\nSaved augmented CSV to:", output_csv_path)
        print("\nSummary statistics:")
        print(summary)

        if args.make_plots:
            if args.output_prefix is None:
                output_directory = output_csv_path.parent
                file_prefix = output_csv_path.stem
            else:
                output_prefix_path = Path(args.output_prefix)
                output_directory = output_prefix_path.parent if output_prefix_path.parent != Path("") else Path(".")
                file_prefix = output_prefix_path.name
            output_directory.mkdir(parents=True, exist_ok=True)
            written_paths = plot_bulkiness(
                output_dataframe=output_dataframe,
                output_dir=output_directory,
                file_prefix=file_prefix,
            )
            print("\nBulkiness plots saved to:", output_directory)
            for written_path in written_paths:
                print(" ", written_path)
        return

    profile_dataframe, summary = compute_nearest_neighbor_distance_profiles(
        points_dataframe=points_dataframe,
        max_neighbors=args.max_neighbors,
    )
    pairwise_csv_path = (
        Path(args.pairwise_output_csv)
        if args.pairwise_output_csv
        else csv_path.with_name(csv_path.stem + "_nearest_neighbor_distance_profiles.csv")
    )
    profile_dataframe.to_csv(pairwise_csv_path, index=False)

    profile_png_path = csv_path.with_name(csv_path.stem + "_nearest_neighbor_distance_profiles.png")
    plot_nearest_neighbor_profiles(
        profile_dataframe=profile_dataframe,
        max_neighbors=int(summary["max_neighbors_used"]),
        output_png_path=profile_png_path,
    )

    print("\nSaved nearest-neighbor profile CSV to:", pairwise_csv_path)
    print("Saved nearest-neighbor profile figure to:", profile_png_path)
    print("\nNearest-neighbor distance summary:")
    print(summary)


if __name__ == "__main__":
    main(sys.argv[1:])
