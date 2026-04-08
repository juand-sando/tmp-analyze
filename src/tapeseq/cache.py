from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np

from .models import AnalysisConfig, SequenceSignal
from .run_logging import log_cache_event


CACHE_VERSION = 1
SCORING_MATRIX_NAME = "BLOSUM45"


def sequence_sha256(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def safe_accession_label(accession: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9._-]+", "_", accession).strip("_")
    return sanitized or "accession"


def cache_entry_dir(cache_dir: Path, accession: str, sequence: str) -> Path:
    sequence_hash = sequence_sha256(sequence)
    return cache_dir / f"{safe_accession_label(accession)}__{sequence_hash[:16]}"


def expected_cache_metadata(accession: str, sequence: str, config: AnalysisConfig) -> dict[str, object]:
    return {
        "cache_version": CACHE_VERSION,
        "accession": accession,
        "sequence_hash": sequence_sha256(sequence),
        "sequence_length": len(sequence),
        "smooth_iterations": config.smooth_iterations,
        "num_shuffles": config.num_shuffles,
        "seed": config.seed,
        "scoring_matrix": SCORING_MATRIX_NAME,
    }


def metadata_matches(actual: dict[str, object], expected: dict[str, object]) -> bool:
    for key, expected_value in expected.items():
        if actual.get(key) != expected_value:
            return False
    return True


def load_cached_full_sequence_signal(
    accession: str,
    full_sequence: str,
    config: AnalysisConfig,
    *,
    need_raw_matrix: bool,
    need_smoothed_matrix: bool,
) -> SequenceSignal | None:
    if config.rebuild_cache:
        log_cache_event("skip", config.cache_dir, detail="rebuild requested")
        return None

    entry_dir = cache_entry_dir(config.cache_dir, accession, full_sequence)
    metadata_path = entry_dir / "metadata.json"
    if not metadata_path.exists():
        log_cache_event("miss", entry_dir, detail="metadata missing")
        return None

    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except Exception:
        log_cache_event("miss", entry_dir, detail="metadata unreadable")
        return None

    expected = expected_cache_metadata(accession, full_sequence, config)
    if not metadata_matches(metadata, expected):
        log_cache_event("miss", entry_dir, detail="metadata mismatch")
        return None

    required_files = [
        entry_dir / "projection_full_raw.npy",
        entry_dir / "projection_full_standardized.npy",
    ]
    if need_raw_matrix:
        required_files.append(entry_dir / "raw_matrix_full.npy")
    if need_smoothed_matrix:
        required_files.append(entry_dir / "smoothed_matrix_full.npy")

    if not all(path.exists() for path in required_files):
        log_cache_event("miss", entry_dir, detail="required artifact missing")
        return None

    signal = SequenceSignal(accession=accession, full_sequence=full_sequence)
    signal.projection_full_raw = np.load(entry_dir / "projection_full_raw.npy")
    signal.projection_full_standardized = np.load(entry_dir / "projection_full_standardized.npy")
    signal.full_baseline_mean = float(metadata["full_baseline_mean"])
    signal.full_baseline_std = float(metadata["full_baseline_std"])

    if need_raw_matrix:
        signal.raw_matrix_full = np.load(entry_dir / "raw_matrix_full.npy").astype(float, copy=False)

    if need_smoothed_matrix:
        signal.smoothed_matrix_full = np.load(entry_dir / "smoothed_matrix_full.npy").astype(float, copy=False)

    log_cache_event("hit", entry_dir)
    return signal


def save_cached_full_sequence_signal(signal: SequenceSignal, config: AnalysisConfig) -> None:
    if (
        signal.raw_matrix_full is None
        or signal.smoothed_matrix_full is None
        or signal.projection_full_raw is None
        or signal.projection_full_standardized is None
        or signal.full_baseline_mean is None
        or signal.full_baseline_std is None
    ):
        return

    entry_dir = cache_entry_dir(config.cache_dir, signal.accession, signal.full_sequence)
    entry_dir.mkdir(parents=True, exist_ok=True)
    overwrite = any((entry_dir / name).exists() for name in [
        "raw_matrix_full.npy",
        "smoothed_matrix_full.npy",
        "projection_full_raw.npy",
        "projection_full_standardized.npy",
        "metadata.json",
    ])

    np.save(entry_dir / "raw_matrix_full.npy", np.rint(signal.raw_matrix_full).astype(np.int16))
    np.save(entry_dir / "smoothed_matrix_full.npy", signal.smoothed_matrix_full.astype(np.float32))
    np.save(entry_dir / "projection_full_raw.npy", signal.projection_full_raw.astype(np.float32))
    np.save(entry_dir / "projection_full_standardized.npy", signal.projection_full_standardized.astype(np.float32))

    metadata = expected_cache_metadata(signal.accession, signal.full_sequence, config)
    metadata["full_baseline_mean"] = float(signal.full_baseline_mean)
    metadata["full_baseline_std"] = float(signal.full_baseline_std)

    with (entry_dir / "metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
    log_cache_event("overwrite" if overwrite else "write", entry_dir)
