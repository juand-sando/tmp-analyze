from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np


@dataclass
class AnalysisConfig:
    mode: str
    input_fasta: Path
    output_dir: Path
    trimmings_csv: Path | None = None
    repeats_fasta: Path | None = None
    smooth_iterations: int = 20
    fft_min: int = 3
    fft_max: int = 19
    min_repeat_length: int = 5
    exact_only: bool = False
    write_svg: bool = False
    num_shuffles: int = 1
    seed: int = 42


@dataclass
class SequenceSignal:
    accession: str
    full_sequence: str
    raw_matrix_full: np.ndarray | None = None
    smoothed_matrix_full: np.ndarray | None = None
    projection_full_raw: np.ndarray | None = None
    projection_full_standardized: np.ndarray | None = None
    trimmed_sequence: str | None = None
    trimming_start_index: int | None = None
    trimming_end_index: int | None = None
    raw_matrix_trimmed: np.ndarray | None = None
    smoothed_matrix_trimmed: np.ndarray | None = None
    projection_trimmed_raw: np.ndarray | None = None
    projection_trimmed_standardized: np.ndarray | None = None
    full_baseline_mean: float | None = None
    full_baseline_std: float | None = None


@dataclass
class PeriodicityResult:
    accession: str
    trimmed_length: int
    best_frequency: int
    best_power: float
    period: float
    repeat_count: int
    repeat_sequences: list[str] = field(default_factory=list)


@dataclass
class SuppliedRepeatScore:
    accession: str
    suffix: str
    repeat_sequence: str
    match_method: str
    window_start: int | None
    window_end: int | None
    window_length: int | None
    score: float | None
