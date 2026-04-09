from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from .analysis import compute_all_pairwise_distances, compute_bulkiness_measures, summarize_pairwise_distances
from .io import load_points
from .plotting import plot_bulkiness, plot_pairwise_distance_histogram


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="Compute bulkiness or pairwise-distance histogram from a 3D C-alpha CSV."
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
        "--radii",
        type=float,
        nargs="+",
        default=[10.0],
        help="One or more radii for neighbor counts.",
    )
    parser.add_argument(
        "--k_neighbors",
        type=int,
        default=10,
        help="k for kNN mean-distance bulkiness.",
    )
    parser.add_argument(
        "--subtract_baseline",
        action="store_true",
        help="Subtract expected neighbors for a straight chain.",
    )
    parser.add_argument(
        "--baseline_spacing",
        type=float,
        default=3.8,
        help="Spacing used for the straight-chain baseline.",
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
        "--histogram_bins",
        type=int,
        default=40,
        help="Number of bins for the pairwise distance histogram.",
    )
    parser.add_argument(
        "--pairwise_output_csv",
        type=str,
        default=None,
        help="Output CSV path for pairwise mode. Default: <input>_pairwise_distances.csv",
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
            radii=args.radii,
            k_neighbors=args.k_neighbors,
            subtract_baseline=args.subtract_baseline,
            baseline_spacing=args.baseline_spacing,
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
                radii=args.radii,
                subtract_baseline=args.subtract_baseline,
                output_dir=output_directory,
                file_prefix=file_prefix,
            )
            print("\nBulkiness plots saved to:", output_directory)
            for written_path in written_paths:
                print(" ", written_path)
        return

    all_pairwise_distances = compute_all_pairwise_distances(points_dataframe)
    pairwise_csv_path = (
        Path(args.pairwise_output_csv)
        if args.pairwise_output_csv
        else csv_path.with_name(csv_path.stem + "_pairwise_distances.csv")
    )
    pd.DataFrame({"pairwise_distance": all_pairwise_distances}).to_csv(pairwise_csv_path, index=False)

    histogram_png_path = csv_path.with_name(csv_path.stem + "_pairwise_distance_hist.png")
    plot_pairwise_distance_histogram(
        all_pairwise_distances=all_pairwise_distances,
        histogram_bins=args.histogram_bins,
        output_png_path=histogram_png_path,
    )

    print("\nSaved pairwise distances CSV to:", pairwise_csv_path)
    print("Saved pairwise distance histogram to:", histogram_png_path)
    print("\nPairwise distance summary:")
    print(summarize_pairwise_distances(all_pairwise_distances, len(points_dataframe)))


if __name__ == "__main__":
    main(sys.argv[1:])
