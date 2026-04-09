from __future__ import annotations

import argparse

from .models import BulkinessConfig


FORMATTER_CLASS = argparse.ArgumentDefaultsHelpFormatter


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-i",
        "--input-file",
        "--input-csv",
        "--input-structure",
        dest="input_path_flag",
        metavar="INPUT_FILE",
        help="Input CSV, PDB, or mmCIF file.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        dest="output_dir_flag",
        metavar="OUTPUT_DIR",
        help="Output directory.",
    )
    parser.add_argument("legacy_input_csv", nargs="?", help=argparse.SUPPRESS)
    parser.add_argument("legacy_output_dir", nargs="?", help=argparse.SUPPRESS)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="tapebulk structure analysis toolkit.",
        formatter_class=FORMATTER_CLASS,
    )
    subparsers = parser.add_subparsers(dest="mode", required=True)

    bulkiness = subparsers.add_parser(
        "bulkiness",
        help="Compute neighbor-count and kNN bulkiness measures from 3D coordinates.",
        formatter_class=FORMATTER_CLASS,
    )
    add_common_arguments(bulkiness)
    bulkiness.add_argument(
        "--radii",
        type=float,
        nargs="+",
        default=[10.0],
        help="One or more radii for neighbor counts.",
    )
    bulkiness.add_argument(
        "--k-neighbors",
        dest="k_neighbors",
        type=int,
        default=10,
        help="k for mean nearest-neighbor distance bulkiness.",
    )
    bulkiness.add_argument(
        "--subtract-baseline",
        action="store_true",
        help="Subtract the expected straight-chain neighbor baseline.",
    )
    bulkiness.add_argument(
        "--baseline-spacing",
        type=float,
        default=3.8,
        help="Residue spacing used for the straight-chain baseline.",
    )
    bulkiness.add_argument(
        "--make-plots",
        action="store_true",
        help="Write bulkiness plots alongside tabular outputs.",
    )

    pairwise = subparsers.add_parser(
        "pairwise-histogram",
        help="Compute all pairwise distances and write a histogram.",
        formatter_class=FORMATTER_CLASS,
    )
    add_common_arguments(pairwise)
    pairwise.add_argument(
        "--histogram-bins",
        type=int,
        default=40,
        help="Number of bins in the pairwise-distance histogram.",
    )

    return parser


def resolve_input_output_args(parser: argparse.ArgumentParser, args) -> tuple[str, str]:
    named_input = getattr(args, "input_path_flag", None)
    named_output = getattr(args, "output_dir_flag", None)
    legacy_input = getattr(args, "legacy_input_csv", None)
    legacy_output = getattr(args, "legacy_output_dir", None)

    if named_input is not None and legacy_input is not None:
        parser.error("provide either the named input flag or the legacy positional input, not both.")
    if named_output is not None and legacy_output is not None:
        parser.error("provide either --output-dir or the legacy positional output, not both.")

    input_value = named_input if named_input is not None else legacy_input
    output_value = named_output if named_output is not None else legacy_output

    if input_value is None or output_value is None:
        parser.error("the following arguments are required: --input-file, --output-dir")

    return input_value, output_value


def parse_args(argv: list[str] | None = None) -> BulkinessConfig:
    parser = build_parser()
    args = parser.parse_args(argv)
    input_path, output_dir = resolve_input_output_args(parser, args)

    config = BulkinessConfig(
        mode=args.mode,
        input_path=input_path,
        output_dir=output_dir,
        radii=list(getattr(args, "radii", [10.0])),
        k_neighbors=getattr(args, "k_neighbors", 10),
        subtract_baseline=getattr(args, "subtract_baseline", False),
        baseline_spacing=getattr(args, "baseline_spacing", 3.8),
        histogram_bins=getattr(args, "histogram_bins", 40),
        make_plots=getattr(args, "make_plots", False),
    )
    validate_config(config)
    return config


def validate_config(config: BulkinessConfig) -> None:
    if config.k_neighbors < 1:
        raise SystemExit("ERROR: --k-neighbors must be >= 1.")
    if any(radius <= 0 for radius in config.radii):
        raise SystemExit("ERROR: all --radii values must be > 0.")
    if config.baseline_spacing <= 0:
        raise SystemExit("ERROR: --baseline-spacing must be > 0.")
    if config.histogram_bins < 1:
        raise SystemExit("ERROR: --histogram-bins must be >= 1.")
