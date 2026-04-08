from __future__ import annotations

import argparse

from .models import AnalysisConfig


FORMATTER_CLASS = argparse.ArgumentDefaultsHelpFormatter


def add_common_signal_arguments(
    parser: argparse.ArgumentParser,
    *,
    include_trimmings: bool = True,
    include_write_svg: bool = True,
) -> None:
    parser.add_argument("input_fasta", help="FASTA of full sequences.")
    parser.add_argument("output_dir", help="Output directory.")
    if include_trimmings:
        parser.add_argument(
            "-t",
            "--trimmings",
            dest="trimmings_csv",
            help="CSV with trimming rules. Rules can be keyed by sequence length or accession.",
        )
    parser.add_argument(
        "--cache-dir",
        default=".tapeseq_cache",
        help="Directory for reusable cached full-sequence analysis artifacts.",
    )
    parser.add_argument(
        "--rebuild-cache",
        action="store_true",
        help="Ignore any matching cached artifacts and recompute them from scratch.",
    )
    parser.add_argument(
        "--smooth-iters",
        dest="smooth_iterations",
        type=int,
        default=20,
        help="Number of smoothing iterations applied to the full self-similarity matrix before projection.",
    )
    if include_write_svg:
        parser.add_argument(
            "--write-svg",
            action="store_true",
            help="Also write SVG versions of the smoothed heatmaps in addition to PNG output.",
        )
    parser.add_argument(
        "--num-shuffles",
        type=int,
        default=1,
        help="Number of shuffled baseline replicates used to standardize the full projection for each accession.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Master random seed used to derive reproducible accession-specific shuffle seeds.",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="tapeseq sequence analysis toolkit.",
        formatter_class=FORMATTER_CLASS,
    )
    subparsers = parser.add_subparsers(dest="mode", required=True)

    exploratory = subparsers.add_parser(
        "exploratory",
        help="Write raw and smoothed heatmaps plus projections for each accession.",
        formatter_class=FORMATTER_CLASS,
    )
    add_common_signal_arguments(exploratory)

    periodicity = subparsers.add_parser(
        "periodicity",
        help="Run Fourier analysis on the focused sequence window and write a periodicity summary.",
        formatter_class=FORMATTER_CLASS,
    )
    add_common_signal_arguments(periodicity)
    periodicity.add_argument(
        "--fft-min",
        type=int,
        default=3,
        help="Minimum FFT frequency index considered when searching for the dominant periodicity peak.",
    )
    periodicity.add_argument(
        "--fft-max",
        type=int,
        default=19,
        help="Maximum FFT frequency index considered when searching for the dominant periodicity peak.",
    )

    repeat_analysis = subparsers.add_parser(
        "repeat-analysis",
        help="Score supplied repeats against the full-sequence projection.",
        formatter_class=FORMATTER_CLASS,
    )
    add_common_signal_arguments(repeat_analysis, include_trimmings=False, include_write_svg=False)
    repeat_analysis.add_argument(
        "-r",
        "--repeats-fasta",
        dest="repeats_fasta",
        help="FASTA of supplied repeats.",
        required=True,
    )
    repeat_analysis.add_argument(
        "--min-repeat-len",
        dest="min_repeat_length",
        type=int,
        default=5,
        help="Minimum supplied repeat length to score.",
    )
    repeat_analysis.add_argument(
        "--exact-only",
        action="store_true",
        help="Disable local-alignment fallback and require exact sequence matches for supplied repeats.",
    )

    map_domain = subparsers.add_parser(
        "map-domain",
        help="Project a reference domain through an alignment and generate accession-based trimming rules.",
        formatter_class=FORMATTER_CLASS,
    )
    map_domain.add_argument("input_fasta", help="FASTA of aligned sequences.")
    map_domain.add_argument("output_dir", help="Output directory.")
    map_domain.add_argument(
        "--reference-accession",
        required=True,
        help="Accession in the alignment used as the coordinate anchor for the reference domain.",
    )
    map_domain.add_argument(
        "--reference-domain",
        dest="reference_domain_fasta",
        required=True,
        help="FASTA containing the reference domain sequence to map through the alignment.",
    )

    return parser


def parse_args(argv: list[str] | None = None) -> AnalysisConfig:
    parser = build_parser()
    args = parser.parse_args(argv)
    config = AnalysisConfig(
        mode=args.mode,
        input_fasta=args.input_fasta,
        output_dir=args.output_dir,
        cache_dir=getattr(args, "cache_dir", ".tapeseq_cache"),
        trimmings_csv=getattr(args, "trimmings_csv", None),
        repeats_fasta=getattr(args, "repeats_fasta", None),
        reference_accession=getattr(args, "reference_accession", None),
        reference_domain_fasta=getattr(args, "reference_domain_fasta", None),
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

    if config.mode == "map-domain":
        if config.reference_accession is None:
            raise SystemExit("ERROR: map-domain requires --reference-accession.")
        if config.reference_domain_fasta is None:
            raise SystemExit("ERROR: map-domain requires --reference-domain.")
        return

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
