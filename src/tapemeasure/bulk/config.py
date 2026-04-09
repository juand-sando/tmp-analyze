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
        help="Compute radius-based neighbor counts and a 3-round chain-aware smoothed profile.",
        formatter_class=FORMATTER_CLASS,
    )
    add_common_arguments(bulkiness)
    bulkiness.add_argument(
        "--radius",
        type=float,
        default=10,
        help="Radius in Angstroms used to count spatial neighbors.",
    )
    bulkiness.add_argument(
        "--make-plots",
        action="store_true",
        help="Write bulkiness plots alongside tabular outputs.",
    )

    pairwise = subparsers.add_parser(
        "pairwise-histogram",
        help="Compute nearest-neighbor distance profiles and write a multi-panel figure.",
        formatter_class=FORMATTER_CLASS,
    )
    add_common_arguments(pairwise)
    pairwise.add_argument(
        "--max-neighbors",
        type=int,
        default=15,
        help="Maximum nearest-neighbor rank to analyze.",
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
        radius=getattr(args, "radius", 10.0),
        max_neighbors=getattr(args, "max_neighbors", 15),
        make_plots=getattr(args, "make_plots", False),
    )
    validate_config(config)
    return config


def validate_config(config: BulkinessConfig) -> None:
    if config.radius <= 0:
        raise SystemExit("ERROR: --radius must be > 0.")
    if config.max_neighbors < 1:
        raise SystemExit("ERROR: --max-neighbors must be >= 1.")
