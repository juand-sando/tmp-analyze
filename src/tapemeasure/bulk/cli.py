from __future__ import annotations

from pathlib import Path

from ..shared.run_logging import finalize_run_logging, setup_run_logging, write_output_inventory
from .config import parse_args
from .workflows import run


def main(argv: list[str] | None = None) -> None:
    config = parse_args(argv)
    config.input_path = Path(config.input_path)
    config.output_dir = Path(config.output_dir)

    config.output_dir.mkdir(parents=True, exist_ok=True)
    setup_run_logging(config)

    summary = "run failed"
    try:
        summary = run(config)
        print(summary)
        return
    except Exception as exc:
        from ..shared.run_logging import get_logger

        get_logger().exception(f"error: run failed with exception: {exc}")
        raise
    finally:
        write_output_inventory(config.output_dir)
        finalize_run_logging(summary or "run completed")
