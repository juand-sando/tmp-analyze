from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from ..shared.run_logging import log_file_write
from .models import BulkinessConfig


def metadata_dir(output_dir: Path) -> Path:
    return output_dir / "run_metadata"


def analysis_dir(output_dir: Path) -> Path:
    return output_dir / "global_outputs" / "bulkiness"


def write_analysis_metadata(config: BulkinessConfig) -> None:
    output_directory = metadata_dir(config.output_dir)
    output_directory.mkdir(parents=True, exist_ok=True)

    config_dict = asdict(config)
    serializable_config = {
        key: str(value) if isinstance(value, Path) else value
        for key, value in config_dict.items()
    }

    metadata_path = output_directory / "analysis_parameters.json"
    log_file_write(metadata_path, category="metadata")
    with metadata_path.open("w", encoding="utf-8") as handle:
        json.dump(serializable_config, handle, indent=2, sort_keys=True)


def write_dataframe_csv(dataframe: pd.DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    log_file_write(output_path)
    dataframe.to_csv(output_path, index=False)


def write_summary_json(summary: dict[str, object], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    log_file_write(output_path, category="metadata")
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
