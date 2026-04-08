from __future__ import annotations

from pathlib import Path

from .config import parse_args
from .run_logging import finalize_run_logging, setup_run_logging, write_output_inventory


def main(argv: list[str] | None = None) -> None:
    config = parse_args(argv)
    config.input_fasta = Path(config.input_fasta)
    config.output_dir = Path(config.output_dir)
    config.cache_dir = Path(config.cache_dir)
    if config.trimmings_csv is not None:
        config.trimmings_csv = Path(config.trimmings_csv)
    if config.repeats_fasta is not None:
        config.repeats_fasta = Path(config.repeats_fasta)
    if config.reference_domain_fasta is not None:
        config.reference_domain_fasta = Path(config.reference_domain_fasta)

    config.output_dir.mkdir(parents=True, exist_ok=True)
    setup_run_logging(config)

    if config.mode == "exploratory":
        from .exploratory import run
    elif config.mode == "periodicity":
        from .periodicity import run
    elif config.mode == "repeat-analysis":
        from .repeat_analysis import run
    elif config.mode == "map-domain":
        from .map_domain import run
    else:
        raise SystemExit(f"ERROR: unsupported mode {config.mode}")

    summary = "run failed"
    try:
        summary = run(config)
        return
    except Exception as exc:
        from .run_logging import get_logger

        get_logger().exception(f"error: run failed with exception: {exc}")
        raise
    finally:
        write_output_inventory(config.output_dir)
        finalize_run_logging(summary or "run completed")
