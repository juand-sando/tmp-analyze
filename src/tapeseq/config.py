from __future__ import annotations

import argparse

from .models import AnalysisConfig


def add_common_signal_arguments(
    parser: argparse.ArgumentParser,
    *,
    include_write_svg: bool = True,
) -> None:
    parser.add_argument("input_fasta", help="FASTA of full sequences.")
    parser.add_argument("output_dir", help="Output directory.")
    parser.add_argument("-t", "--trimmings", dest="trimmings_csv", help="CSV with trimming rules.")
    parser.add_argument("--cache-dir", default=".tapeseq_cache", help="Directory for reusable cached analysis artifacts.")
    parser.add_argument("--rebuild-cache", action="store_true", help="Ignore cached artifacts and rebuild them.")
    parser.add_argument("--smooth-iters", dest="smooth_iterations", type=int, default=20)
    if include_write_svg:
        parser.add_argument("--write-svg", action="store_true")
    parser.add_argument("--num-shuffles", type=int, default=1)
    parser.add_argument("--seed", type=int, default=42)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="tapeseq sequence analysis toolkit.")
    subparsers = parser.add_subparsers(dest="mode", required=True)

    exploratory = subparsers.add_parser(
        "exploratory",
        help="Write raw and smoothed heatmaps plus projections for each accession.",
    )
    add_common_signal_arguments(exploratory)

    periodicity = subparsers.add_parser(
        "periodicity",
        help="Run Fourier analysis on the focused sequence window and write a periodicity summary.",
    )
    add_common_signal_arguments(periodicity)
    periodicity.add_argument("--fft-min", type=int, default=3)
    periodicity.add_argument("--fft-max", type=int, default=19)

    repeat_analysis = subparsers.add_parser(
        "repeat-analysis",
        help="Score supplied repeats against the full-sequence projection.",
    )
    add_common_signal_arguments(repeat_analysis, include_write_svg=False)
    repeat_analysis.add_argument(
        "-r",
        "--repeats-fasta",
        dest="repeats_fasta",
        help="FASTA of supplied repeats.",
        required=True,
    )
    repeat_analysis.add_argument("--min-repeat-len", dest="min_repeat_length", type=int, default=5)
    repeat_analysis.add_argument("--exact-only", action="store_true")

    return parser


def parse_args(argv: list[str] | None = None) -> AnalysisConfig:
    parser = build_parser()
    args = parser.parse_args(argv)
    config = AnalysisConfig(
        mode=args.mode,
        input_fasta=args.input_fasta,
        output_dir=args.output_dir,
        cache_dir=args.cache_dir,
        trimmings_csv=getattr(args, "trimmings_csv", None),
        repeats_fasta=getattr(args, "repeats_fasta", None),
        smooth_iterations=getattr(args, "smooth_iterations", 20),
        fft_min=getattr(args, "fft_min", 3),
        fft_max=getattr(args, "fft_max", 19),
        min_repeat_length=getattr(args, "min_repeat_length", 5),
        exact_only=getattr(args, "exact_only", False),
        write_svg=getattr(args, "write_svg", False),
        num_shuffles=getattr(args, "num_shuffles", 1),
        seed=getattr(args, "seed", 42),
        rebuild_cache=getattr(args, "rebuild_cache", False),
    )
    validate_config(config)
    return config


def validate_config(config: AnalysisConfig) -> None:
    if config.mode == "repeat-analysis" and config.repeats_fasta is None:
        raise SystemExit("ERROR: repeat-analysis requires --repeats-fasta.")

    if config.num_shuffles < 1:
        raise SystemExit("ERROR: --num-shuffles must be >= 1.")

    if config.smooth_iterations < 0:
        raise SystemExit("ERROR: --smooth-iters must be >= 0.")

    if config.fft_min < 1 or config.fft_max < 1:
        raise SystemExit("ERROR: --fft-min and --fft-max must be >= 1.")

    if config.fft_max < config.fft_min:
        raise SystemExit("ERROR: --fft-max must be >= --fft-min.")

    if config.min_repeat_length < 1:
        raise SystemExit("ERROR: --min-repeat-len must be >= 1.")
