from __future__ import annotations

import json
import logging
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path


LOGGER_NAME = "tapeseq.run"


def get_logger() -> logging.Logger:
    return logging.getLogger(LOGGER_NAME)


def setup_run_logging(config) -> Path:
    metadata_directory = config.output_dir / "run_metadata"
    metadata_directory.mkdir(parents=True, exist_ok=True)
    log_path = metadata_directory / "run.log"

    logger = get_logger()
    logger.setLevel(logging.INFO)
    logger.propagate = False

    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(file_handler)

    run_timestamp = datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")
    logger.info("")
    logger.info("=== tapeseq run start ===")
    logger.info(f"time: {run_timestamp}")
    logger.info(f"mode: {config.mode}")
    logger.info(f"output_dir: {config.output_dir}")

    config_dict = asdict(config)
    serializable_config = {
        key: str(value) if isinstance(value, Path) else value
        for key, value in config_dict.items()
    }
    logger.info("parameters:")
    logger.info(json.dumps(serializable_config, indent=2, sort_keys=True))
    return log_path


def finalize_run_logging(message: str) -> None:
    logger = get_logger()
    logger.info(f"summary: {message}")
    logger.info("=== tapeseq run end ===")


def log_file_write(path: Path, *, category: str = "output") -> None:
    logger = get_logger()
    action = "overwrite" if path.exists() else "write"
    logger.info(f"{action} {category}: {path}")


def log_cache_event(event: str, path: Path, *, detail: str | None = None) -> None:
    logger = get_logger()
    suffix = f" ({detail})" if detail else ""
    logger.info(f"cache {event}: {path}{suffix}")


def write_output_inventory(output_dir: Path) -> Path:
    metadata_directory = output_dir / "run_metadata"
    metadata_directory.mkdir(parents=True, exist_ok=True)
    inventory_path = metadata_directory / "current_state.json"

    inventory_records: list[dict[str, object]] = []
    for path in sorted(p for p in output_dir.rglob("*") if p.is_file()):
        relative_path = path.relative_to(output_dir)
        inventory_records.append(
            {
                "path": str(relative_path),
                "size_bytes": path.stat().st_size,
                "modified_time": datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).astimezone().isoformat(timespec="seconds"),
            }
        )

    log_file_write(inventory_path, category="metadata")
    with inventory_path.open("w", encoding="utf-8") as handle:
        json.dump({"files": inventory_records}, handle, indent=2)

    logger = get_logger()
    logger.info(f"current_state_file_count: {len(inventory_records)}")
    logger.info(f"current_state_path: {inventory_path}")
    return inventory_path
