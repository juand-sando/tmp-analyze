from __future__ import annotations

from pathlib import Path

from .config import parse_args


def main(argv: list[str] | None = None) -> None:
    config = parse_args(argv)
    config.input_fasta = Path(config.input_fasta)
    config.output_dir = Path(config.output_dir)
    if config.trimmings_csv is not None:
        config.trimmings_csv = Path(config.trimmings_csv)
    if config.repeats_fasta is not None:
        config.repeats_fasta = Path(config.repeats_fasta)

    if config.mode == "exploratory":
        from .exploratory import run
    elif config.mode == "periodicity":
        from .periodicity import run
    elif config.mode == "repeat-analysis":
        from .repeat_analysis import run
    else:
        raise SystemExit(f"ERROR: unsupported mode {config.mode}")

    run(config)
