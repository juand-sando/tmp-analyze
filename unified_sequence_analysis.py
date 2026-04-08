#!/usr/bin/env python3
"""
Unified sequence analysis pipeline.

Supported analysis branches:
- full-sequence heatmaps
- scrambled full-sequence heatmaps
- trimmed repeat discovery from projection
- scoring supplied repeats against the full-sequence projection

Output structure:
output_dir/
├─ per_accession/
│  └─ ACCESSION/
│     ├─ full_signal/
│     ├─ trimmed_discovery/
│     └─ supplied_repeat_scoring/
├─ global_outputs/
│  ├─ trimmed_discovery/
│  │  ├─ trimmed_sequences.fasta
│  │  └─ discovery_results_summary.csv
│  └─ supplied_repeat_scoring/
│     └─ supplied_repeat_scores_summary.csv
└─ run_metadata/
   └─ analysis_parameters.json

Notes:
- Full-sequence baselines are always computed from the FULL sequence, even when the
  downstream analysis acts on the trimmed signal.
- No discovery_repeat_window_values.csv is produced.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO, pairwise2
from Bio.Align import substitution_matrices as substitution_matrices

mpl.rcParams["svg.fonttype"] = "none"

BLOSUM45 = substitution_matrices.load("BLOSUM45")
REPEAT_HEADER_PATTERN = re.compile(r"^(?P<accession>.+?)_R(?P<suffix>\d+|F)$")


# -----------------------------------------------------------------------------
# Data classes
# -----------------------------------------------------------------------------

@dataclass
class AnalysisConfig:
    input_fasta: Path
    output_dir: Path
    trimmings_csv: Path | None = None
    repeats_fasta: Path | None = None

    full_heatmaps: bool = False
    scrambled_heatmaps: bool = False
    repeat_discovery: bool = False
    score_supplied_repeats: bool = False

    signal_source: str = "raw"  # raw | standardized
    baseline_mode: str = "none"  # none | global_projection
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
class RepeatDiscoveryResult:
    accession: str
    trimmed_length: int
    best_frequency: int
    best_power: float
    period: float
    repeat_count: int
    repeat_sequences: List[str] = field(default_factory=list)


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


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------


def parse_args() -> AnalysisConfig:
    parser = argparse.ArgumentParser(description="Unified sequence analysis pipeline.")

    parser.add_argument("input_fasta", type=Path, help="FASTA of full sequences.")
    parser.add_argument("output_dir", type=Path, help="Output directory.")

    parser.add_argument("-t", "--trimmings", dest="trimmings_csv", type=Path, help="CSV with trimming rules.")
    parser.add_argument("-r", "--repeats-fasta", dest="repeats_fasta", type=Path, help="FASTA of supplied repeats.")

    parser.add_argument("--full-heatmaps", action="store_true", help="Write full-sequence heatmaps.")
    parser.add_argument("--scrambled-heatmaps", action="store_true", help="Write scrambled-control heatmaps.")
    parser.add_argument("--repeat-discovery", action="store_true", help="Run trimmed repeat discovery.")
    parser.add_argument("--score-supplied-repeats", action="store_true", help="Score supplied repeats on the full projection.")

    parser.add_argument("--signal-source", choices=["raw", "standardized"], default="raw")
    parser.add_argument("--baseline-mode", choices=["none", "global_projection"], default="none")

    parser.add_argument("--smooth-iters", dest="smooth_iterations", type=int, default=20)
    parser.add_argument("--fft-min", type=int, default=3)
    parser.add_argument("--fft-max", type=int, default=19)
    parser.add_argument("--min-repeat-len", dest="min_repeat_length", type=int, default=5)
    parser.add_argument("--exact-only", action="store_true")
    parser.add_argument("--write-svg", action="store_true")
    parser.add_argument("--num-shuffles", type=int, default=1)
    parser.add_argument("--seed", type=int, default=42)

    args = parser.parse_args()

    return AnalysisConfig(
        input_fasta=args.input_fasta,
        output_dir=args.output_dir,
        trimmings_csv=args.trimmings_csv,
        repeats_fasta=args.repeats_fasta,
        full_heatmaps=args.full_heatmaps,
        scrambled_heatmaps=args.scrambled_heatmaps,
        repeat_discovery=args.repeat_discovery,
        score_supplied_repeats=args.score_supplied_repeats,
        signal_source=args.signal_source,
        baseline_mode=args.baseline_mode,
        smooth_iterations=args.smooth_iterations,
        fft_min=args.fft_min,
        fft_max=args.fft_max,
        min_repeat_length=args.min_repeat_length,
        exact_only=args.exact_only,
        write_svg=args.write_svg,
        num_shuffles=args.num_shuffles,
        seed=args.seed,
    )


def validate_config(config: AnalysisConfig) -> None:
    if not any([
        config.full_heatmaps,
        config.scrambled_heatmaps,
        config.repeat_discovery,
        config.score_supplied_repeats,
    ]):
        config.full_heatmaps = True

    if config.repeat_discovery and config.trimmings_csv is None:
        raise SystemExit("ERROR: --repeat-discovery requires --trimmings.")

    if config.score_supplied_repeats and config.repeats_fasta is None:
        raise SystemExit("ERROR: --score-supplied-repeats requires --repeats-fasta.")

    if config.signal_source == "standardized" and config.baseline_mode == "none":
        raise SystemExit("ERROR: --signal-source standardized requires --baseline-mode global_projection.")

    if config.num_shuffles < 1:
        raise SystemExit("ERROR: --num-shuffles must be >= 1.")

    if config.fft_min < 1 or config.fft_max < 1:
        raise SystemExit("ERROR: --fft-min and --fft-max must be >= 1.")


# -----------------------------------------------------------------------------
# Paths and metadata
# -----------------------------------------------------------------------------


def per_accession_base_dir(output_dir: Path, accession: str) -> Path:
    return output_dir / "per_accession" / accession


def full_signal_dir(output_dir: Path, accession: str) -> Path:
    return per_accession_base_dir(output_dir, accession) / "full_signal"


def trimmed_discovery_dir(output_dir: Path, accession: str) -> Path:
    return per_accession_base_dir(output_dir, accession) / "trimmed_discovery"


def supplied_repeat_scoring_dir(output_dir: Path, accession: str) -> Path:
    return per_accession_base_dir(output_dir, accession) / "supplied_repeat_scoring"


def global_trimmed_discovery_dir(output_dir: Path) -> Path:
    return output_dir / "global_outputs" / "trimmed_discovery"


def global_supplied_repeat_scoring_dir(output_dir: Path) -> Path:
    return output_dir / "global_outputs" / "supplied_repeat_scoring"


def metadata_dir(output_dir: Path) -> Path:
    return output_dir / "run_metadata"


def write_analysis_metadata(config: AnalysisConfig) -> None:
    metadata_directory = metadata_dir(config.output_dir)
    metadata_directory.mkdir(parents=True, exist_ok=True)

    config_dict = asdict(config)
    serializable_config = {
        key: str(value) if isinstance(value, Path) else value
        for key, value in config_dict.items()
    }

    with (metadata_directory / "analysis_parameters.json").open("w", encoding="utf-8") as handle:
        json.dump(serializable_config, handle, indent=2)


# -----------------------------------------------------------------------------
# FASTA helpers
# -----------------------------------------------------------------------------


def extract_accession(record) -> str:
    return str(record.description or record.id).split()[0]



def load_full_sequences(fasta_path: Path) -> Dict[str, str]:
    sequences_by_accession: Dict[str, str] = {}
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        accession = extract_accession(record)
        sequences_by_accession[accession] = str(record.seq).upper()
    return sequences_by_accession



def parse_repeat_header(header: str) -> Tuple[str | None, str | None]:
    token = str(header).split()[0]
    match = REPEAT_HEADER_PATTERN.match(token)
    if not match:
        return None, None
    return match.group("accession"), match.group("suffix")



def load_supplied_repeats(repeats_fasta: Path) -> List[Tuple[str, str, str]]:
    repeat_records: List[Tuple[str, str, str]] = []
    for record in SeqIO.parse(str(repeats_fasta), "fasta"):
        accession, suffix = parse_repeat_header(str(record.description or record.id))
        if accession is None or suffix is None:
            continue
        repeat_sequence = re.sub(r"[\s\-]", "", str(record.seq).upper())
        repeat_records.append((accession, suffix, repeat_sequence))
    return repeat_records


# -----------------------------------------------------------------------------
# Trimming
# -----------------------------------------------------------------------------


def load_trimmings(csv_path: Path) -> Dict[int, Tuple[int, int]]:
    trimming_map: Dict[int, Tuple[int, int]] = {}
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise SystemExit(f"ERROR: {csv_path} has no header row.")
        normalized_headers = [header.lstrip("\ufeff").strip() for header in reader.fieldnames]
        reader.fieldnames = normalized_headers

        required_columns = {"tmp_length", "start", "end"}
        if not required_columns.issubset(set(normalized_headers)):
            raise SystemExit(f"ERROR: {csv_path} must contain columns {sorted(required_columns)}.")

        for row in reader:
            try:
                length_value = int((row["tmp_length"] or "").strip())
                start_value = int((row["start"] or "").strip())
                end_value = int((row["end"] or "").strip())
            except Exception:
                continue
            if start_value >= 1:
                trimming_map[length_value] = (start_value, end_value)

    if not trimming_map:
        raise SystemExit(f"ERROR: No valid trimming rows found in {csv_path}.")
    return trimming_map



def trimming_bounds(full_length: int, trim_map: Dict[int, Tuple[int, int]]) -> Tuple[int, int] | None:
    if full_length not in trim_map:
        return None

    start_1_based, end_specification = trim_map[full_length]
    start_index = start_1_based - 1

    if end_specification <= 0:
        end_index = full_length + end_specification
    else:
        end_index = end_specification

    if not (0 <= start_index < end_index <= full_length):
        return None

    return start_index, end_index


# -----------------------------------------------------------------------------
# Core scoring and signal generation
# -----------------------------------------------------------------------------


def score_pair(residue_a: str, residue_b: str, default: int = 0) -> int:
    residue_a = residue_a.upper()
    residue_b = residue_b.upper()
    try:
        return BLOSUM45[(residue_a, residue_b)]
    except KeyError:
        try:
            return BLOSUM45[(residue_b, residue_a)]
        except KeyError:
            return default



def build_self_matrix(sequence: str) -> np.ndarray:
    sequence_length = len(sequence)
    matrix = np.zeros((sequence_length, sequence_length), dtype=float)
    for row_index, residue_a in enumerate(sequence):
        for column_index in range(row_index, sequence_length):
            score_value = score_pair(residue_a, sequence[column_index])
            matrix[row_index, column_index] = score_value
            if column_index != row_index:
                matrix[column_index, row_index] = score_value
    return matrix



def smooth_matrix(matrix: np.ndarray, iterations: int = 20) -> np.ndarray:
    row_count, column_count = matrix.shape
    for _ in range(iterations):
        new_matrix = np.copy(matrix)
        for row_index in range(1, row_count):
            for column_index in range(1, column_count):
                new_matrix[row_index, column_index] = 0.5 * matrix[row_index, column_index] + 0.5 * max(
                    matrix[row_index - 1, column_index - 1],
                    matrix[row_index - 1, column_index] / 2,
                    matrix[row_index, column_index - 1] / 2,
                )
        matrix = np.copy(new_matrix)
        for row_index in range(row_count - 2, -1, -1):
            for column_index in range(column_count - 2, -1, -1):
                matrix[row_index, column_index] = 0.5 * new_matrix[row_index, column_index] + 0.5 * max(
                    new_matrix[row_index + 1, column_index + 1],
                    new_matrix[row_index + 1, column_index] / 2,
                    new_matrix[row_index, column_index + 1] / 2,
                )
    return matrix



def project_rowwise_max(matrix: np.ndarray) -> np.ndarray:
    masked_matrix = matrix.copy()
    np.fill_diagonal(masked_matrix, -np.inf)
    return np.max(masked_matrix, axis=1)



def scrambled_sequence(sequence: str, seed: int) -> str:
    rng = random.Random(seed)
    residues = list(sequence)
    rng.shuffle(residues)
    return "".join(residues)



def build_projection_baseline(
    sequence: str,
    smooth_iterations: int,
    num_shuffles: int,
    seed: int,
) -> Tuple[float, float]:
    randomized_projection_values: List[np.ndarray] = []

    for shuffle_index in range(num_shuffles):
        shuffled_sequence = scrambled_sequence(sequence, seed + shuffle_index)
        randomized_matrix = build_self_matrix(shuffled_sequence)
        randomized_smoothed_matrix = smooth_matrix(randomized_matrix, iterations=smooth_iterations)
        randomized_projection = project_rowwise_max(randomized_smoothed_matrix)
        randomized_projection_values.append(randomized_projection)

    pooled_values = np.concatenate(randomized_projection_values)
    baseline_mean = float(np.mean(pooled_values))
    baseline_std = float(np.std(pooled_values))
    return baseline_mean, baseline_std



def standardize_signal(signal: np.ndarray, baseline_mean: float, baseline_std: float) -> np.ndarray:
    safe_standard_deviation = baseline_std if baseline_std > 1e-12 else 1.0
    return (signal - baseline_mean) / safe_standard_deviation



def build_full_sequence_signal(accession: str, full_sequence: str, config: AnalysisConfig) -> SequenceSignal:
    signal = SequenceSignal(accession=accession, full_sequence=full_sequence)

    if len(full_sequence) < 2:
        return signal

    signal.raw_matrix_full = build_self_matrix(full_sequence)
    signal.smoothed_matrix_full = smooth_matrix(signal.raw_matrix_full.copy(), iterations=config.smooth_iterations)
    signal.projection_full_raw = project_rowwise_max(signal.smoothed_matrix_full)

    if config.baseline_mode == "global_projection":
        baseline_mean, baseline_std = build_projection_baseline(
            sequence=full_sequence,
            smooth_iterations=config.smooth_iterations,
            num_shuffles=config.num_shuffles,
            seed=config.seed,
        )
        signal.full_baseline_mean = baseline_mean
        signal.full_baseline_std = baseline_std
        signal.projection_full_standardized = standardize_signal(
            signal.projection_full_raw,
            baseline_mean,
            baseline_std,
        )

    return signal



def add_trimmed_signal(signal: SequenceSignal, trim_map: Dict[int, Tuple[int, int]], config: AnalysisConfig) -> SequenceSignal:
    bounds = trimming_bounds(len(signal.full_sequence), trim_map)
    if bounds is None:
        return signal

    start_index, end_index = bounds
    trimmed_sequence = signal.full_sequence[start_index:end_index]
    if len(trimmed_sequence) < 2:
        return signal

    signal.trimmed_sequence = trimmed_sequence
    signal.trimming_start_index = start_index
    signal.trimming_end_index = end_index
    signal.raw_matrix_trimmed = build_self_matrix(trimmed_sequence)
    signal.smoothed_matrix_trimmed = smooth_matrix(signal.raw_matrix_trimmed.copy(), iterations=config.smooth_iterations)
    signal.projection_trimmed_raw = project_rowwise_max(signal.smoothed_matrix_trimmed)

    # User requested: baseline is always computed from the FULL sequence.
    if config.baseline_mode == "global_projection" and signal.full_baseline_mean is not None and signal.full_baseline_std is not None:
        signal.projection_trimmed_standardized = standardize_signal(
            signal.projection_trimmed_raw,
            signal.full_baseline_mean,
            signal.full_baseline_std,
        )

    return signal



def choose_projection(signal: SequenceSignal, use_trimmed: bool, signal_source: str) -> np.ndarray | None:
    if use_trimmed:
        if signal_source == "standardized":
            return signal.projection_trimmed_standardized
        return signal.projection_trimmed_raw

    if signal_source == "standardized":
        return signal.projection_full_standardized
    return signal.projection_full_raw


# -----------------------------------------------------------------------------
# Plotting and file writing
# -----------------------------------------------------------------------------


def save_heatmap(matrix: np.ndarray, title: str, out_path: Path) -> None:
    plt.figure()
    plt.imshow(matrix, aspect="auto")
    plt.title(title)
    plt.xlabel("Position")
    plt.ylabel("Position")
    plt.colorbar(label="Score")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()



def save_projection_plot(signal: np.ndarray, title: str, out_path: Path) -> None:
    plt.figure()
    plt.plot(signal)
    plt.title(title)
    plt.xlabel("Position")
    plt.ylabel("Signal")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()



def save_fft_plot(
    power_spectrum: np.ndarray,
    title: str,
    out_path: Path,
    minimum_frequency: int,
    maximum_frequency: int,
    highlight_frequency: int | None = None,
) -> None:
    lower_bound = max(1, minimum_frequency)
    upper_bound = min(maximum_frequency, len(power_spectrum) - 1)
    trimmed_power = power_spectrum[lower_bound:upper_bound + 1]
    frequencies = np.arange(lower_bound, lower_bound + len(trimmed_power))

    plt.figure()
    plt.plot(frequencies, trimmed_power)
    if highlight_frequency is not None:
        plt.axvline(highlight_frequency, linestyle="--")
    plt.title(title)
    plt.xlabel("Frequency index")
    plt.ylabel("Power")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()



def write_wrapped_fasta_record(handle, header: str, sequence: str, line_width: int = 60) -> None:
    handle.write(f">{header}\n")
    for start_index in range(0, len(sequence), line_width):
        handle.write(sequence[start_index:start_index + line_width] + "\n")



def write_matrix_csv(matrix: np.ndarray, sequence: str, out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([""] + list(sequence))
        for residue, row in zip(sequence, matrix):
            writer.writerow([residue] + [f"{float(value):.4f}" for value in row])



def write_full_heatmap_outputs(signal: SequenceSignal, config: AnalysisConfig) -> None:
    if signal.raw_matrix_full is None or signal.smoothed_matrix_full is None:
        return

    output_directory = full_signal_dir(config.output_dir, signal.accession)
    output_directory.mkdir(parents=True, exist_ok=True)

    save_heatmap(
        signal.raw_matrix_full,
        f"{signal.accession} raw matrix (full)",
        output_directory / f"{signal.accession}_full_raw.png",
    )
    save_heatmap(
        signal.smoothed_matrix_full,
        f"{signal.accession} smoothed matrix (full)",
        output_directory / f"{signal.accession}_full_smooth.png",
    )

    if config.write_svg:
        save_heatmap(
            signal.smoothed_matrix_full,
            f"{signal.accession} smoothed matrix (full)",
            output_directory / f"{signal.accession}_full_smooth.svg",
        )



def write_scrambled_heatmap_outputs(signal: SequenceSignal, config: AnalysisConfig) -> None:
    if len(signal.full_sequence) < 2:
        return

    output_directory = full_signal_dir(config.output_dir, signal.accession)
    output_directory.mkdir(parents=True, exist_ok=True)

    shuffled_sequence = scrambled_sequence(signal.full_sequence, config.seed)
    shuffled_raw_matrix = build_self_matrix(shuffled_sequence)
    shuffled_smoothed_matrix = smooth_matrix(shuffled_raw_matrix.copy(), iterations=config.smooth_iterations)

    save_heatmap(
        shuffled_raw_matrix,
        f"{signal.accession} scrambled raw matrix (full)",
        output_directory / f"{signal.accession}_scrambled_raw.png",
    )
    save_heatmap(
        shuffled_smoothed_matrix,
        f"{signal.accession} scrambled smoothed matrix (full)",
        output_directory / f"{signal.accession}_scrambled_smooth.png",
    )


# -----------------------------------------------------------------------------
# Repeat discovery
# -----------------------------------------------------------------------------


def fft_power_and_peak(signal: np.ndarray, minimum_frequency: int, maximum_frequency: int) -> Tuple[np.ndarray, int, float]:
    signal_float = signal.astype(float)
    fft_spectrum = np.fft.fft(signal_float)
    power_spectrum = fft_spectrum.real ** 2 + fft_spectrum.imag ** 2

    lower_bound = max(1, minimum_frequency)
    upper_bound = min(maximum_frequency, len(power_spectrum) - 1)
    if upper_bound < lower_bound:
        return power_spectrum, 0, 0.0

    candidate_indices = np.arange(lower_bound, upper_bound + 1, dtype=int)
    best_index = int(candidate_indices[np.argmax(power_spectrum[candidate_indices])])
    return power_spectrum, best_index, float(power_spectrum[best_index])



def period_from_frequency(signal_length: int, best_frequency: int) -> float:
    if best_frequency <= 0:
        return float("inf")
    return float(signal_length) / float(best_frequency)



def slice_repeats(sequence: str, period: float) -> List[str]:
    if not sequence or not np.isfinite(period) or period < 1.0:
        return []

    window_size = max(1, int(math.floor(period)))
    repeat_sequences: List[str] = []
    start_index = 0

    while start_index < len(sequence):
        remaining_length = len(sequence) - start_index
        if remaining_length < 2 * period:
            end_index = len(sequence)
        else:
            end_index = min(len(sequence), start_index + window_size)

        if end_index <= start_index:
            break

        repeat_sequences.append(sequence[start_index:end_index])

        if end_index == len(sequence):
            break
        start_index += window_size

    return repeat_sequences



def run_repeat_discovery(signal: SequenceSignal, config: AnalysisConfig) -> RepeatDiscoveryResult | None:
    if signal.trimmed_sequence is None:
        return None

    projection = choose_projection(signal, use_trimmed=True, signal_source=config.signal_source)
    if projection is None or projection.size == 0:
        return None

    power_spectrum, best_frequency, best_power = fft_power_and_peak(
        projection,
        minimum_frequency=config.fft_min,
        maximum_frequency=config.fft_max,
    )
    period = period_from_frequency(len(projection), best_frequency)
    repeat_sequences = slice_repeats(signal.trimmed_sequence, period)

    return RepeatDiscoveryResult(
        accession=signal.accession,
        trimmed_length=len(signal.trimmed_sequence),
        best_frequency=best_frequency,
        best_power=best_power,
        period=period,
        repeat_count=len(repeat_sequences),
        repeat_sequences=repeat_sequences,
    )



def write_repeat_discovery_outputs(signal: SequenceSignal, result: RepeatDiscoveryResult, config: AnalysisConfig) -> None:
    output_directory = trimmed_discovery_dir(config.output_dir, signal.accession)
    output_directory.mkdir(parents=True, exist_ok=True)

    if signal.trimmed_sequence is not None:
        with (output_directory / f"{signal.accession}_trimmed.fasta").open("w", encoding="utf-8") as handle:
            write_wrapped_fasta_record(handle, signal.accession, signal.trimmed_sequence)

    if signal.raw_matrix_trimmed is not None:
        save_heatmap(
            signal.raw_matrix_trimmed,
            f"{signal.accession} raw matrix (trimmed)",
            output_directory / f"{signal.accession}_trimmed_raw.png",
        )

    if signal.smoothed_matrix_trimmed is not None:
        save_heatmap(
            signal.smoothed_matrix_trimmed,
            f"{signal.accession} smoothed matrix (trimmed)",
            output_directory / f"{signal.accession}_trimmed_smooth.png",
        )
        if config.write_svg:
            save_heatmap(
                signal.smoothed_matrix_trimmed,
                f"{signal.accession} smoothed matrix (trimmed)",
                output_directory / f"{signal.accession}_trimmed_smooth.svg",
            )

    if signal.smoothed_matrix_trimmed is not None and signal.trimmed_sequence is not None:
        write_matrix_csv(
            signal.smoothed_matrix_trimmed,
            signal.trimmed_sequence,
            output_directory / f"{signal.accession}_trimmed_matrix.csv",
        )

    projection = choose_projection(signal, use_trimmed=True, signal_source=config.signal_source)
    if projection is not None:
        projection_label = "projection"
        save_projection_plot(
            projection,
            f"{signal.accession} {projection_label}",
            output_directory / f"{signal.accession}_projection.png",
        )
        power_spectrum, _, _ = fft_power_and_peak(
            projection,
            minimum_frequency=config.fft_min,
            maximum_frequency=config.fft_max,
        )
        save_fft_plot(
            power_spectrum,
            f"{signal.accession} projection FFT power",
            output_directory / f"{signal.accession}_projection_fft.png",
            minimum_frequency=config.fft_min,
            maximum_frequency=config.fft_max,
            highlight_frequency=result.best_frequency,
        )

    with (output_directory / f"{signal.accession}_repeats_discovered.fasta").open("w", encoding="utf-8") as handle:
        for repeat_index, repeat_sequence in enumerate(result.repeat_sequences, start=1):
            suffix = f"R{repeat_index}" if repeat_index < len(result.repeat_sequences) else "RF"
            write_wrapped_fasta_record(handle, f"{signal.accession}_{suffix}", repeat_sequence)



def write_trimmed_sequences_fasta(signals_by_accession: Dict[str, SequenceSignal], output_dir: Path) -> None:
    output_directory = global_trimmed_discovery_dir(output_dir)
    output_directory.mkdir(parents=True, exist_ok=True)

    trimmed_signals = [signal for signal in signals_by_accession.values() if signal.trimmed_sequence]
    if not trimmed_signals:
        return

    with (output_directory / "trimmed_sequences.fasta").open("w", encoding="utf-8") as handle:
        for signal in trimmed_signals:
            write_wrapped_fasta_record(handle, signal.accession, signal.trimmed_sequence or "")



def write_discovery_results_summary(results: List[RepeatDiscoveryResult], signals_by_accession: Dict[str, SequenceSignal], config: AnalysisConfig) -> None:
    if not results:
        return

    output_directory = global_trimmed_discovery_dir(config.output_dir)
    output_directory.mkdir(parents=True, exist_ok=True)

    summary_path = output_directory / "discovery_results_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        headers = ["accession", "trimmed_length", "max_frequency", "max_power", "period", "repeat_count"]
        if config.baseline_mode == "global_projection":
            headers.extend(["full_baseline_mean", "full_baseline_std"])
        writer.writerow(headers)

        for result in results:
            signal = signals_by_accession[result.accession]
            row = [
                result.accession,
                result.trimmed_length,
                result.best_frequency,
                f"{result.best_power:.6f}",
                f"{result.period:.6f}",
                result.repeat_count,
            ]
            if config.baseline_mode == "global_projection":
                row.extend([
                    "" if signal.full_baseline_mean is None else f"{signal.full_baseline_mean:.6f}",
                    "" if signal.full_baseline_std is None else f"{signal.full_baseline_std:.6f}",
                ])
            writer.writerow(row)


# -----------------------------------------------------------------------------
# Supplied repeat scoring
# -----------------------------------------------------------------------------


def find_repeat_range_exact(full_sequence: str, repeat_sequence: str) -> Tuple[int, int] | None:
    start_index = full_sequence.find(repeat_sequence)
    if start_index < 0:
        return None
    return start_index, start_index + len(repeat_sequence)



def find_repeat_range_align(full_sequence: str, repeat_sequence: str) -> Tuple[int, int] | None:
    alignments = pairwise2.align.localms(full_sequence, repeat_sequence, 2, -1, -5, -1)
    if not alignments:
        return None

    best_alignment = alignments[0]
    if hasattr(best_alignment, "seqA"):
        aligned_full_sequence = best_alignment.seqA
        aligned_repeat_sequence = best_alignment.seqB
        full_start = int(best_alignment.start)
    else:
        aligned_full_sequence, aligned_repeat_sequence, _score, full_start, _end = best_alignment

    full_position = full_start
    mapped_start = None
    mapped_end = None

    for full_character, repeat_character in zip(aligned_full_sequence, aligned_repeat_sequence):
        current_full_position = None
        if full_character != "-":
            current_full_position = full_position
            full_position += 1

        if repeat_character != "-" and current_full_position is not None:
            if mapped_start is None:
                mapped_start = current_full_position
            mapped_end = current_full_position

    if mapped_start is None or mapped_end is None:
        return None

    return mapped_start, mapped_end + 1



def locate_repeat_window(full_sequence: str, repeat_sequence: str, exact_only: bool = False) -> Tuple[Tuple[int, int] | None, str]:
    exact_match = find_repeat_range_exact(full_sequence, repeat_sequence)
    if exact_match is not None:
        return exact_match, "exact"
    if exact_only:
        return None, "not_found"
    alignment_match = find_repeat_range_align(full_sequence, repeat_sequence)
    if alignment_match is not None:
        return alignment_match, "alignment"
    return None, "not_found"



def window_max(signal_vector: np.ndarray, start_index: int, end_index: int) -> float:
    if start_index < 0 or end_index > len(signal_vector) or start_index >= end_index:
        return np.nan
    return float(np.nanmax(signal_vector[start_index:end_index]))



def score_supplied_repeats(
    signals_by_accession: Dict[str, SequenceSignal],
    repeat_records: List[Tuple[str, str, str]],
    config: AnalysisConfig,
) -> List[SuppliedRepeatScore]:
    scored_repeats: List[SuppliedRepeatScore] = []

    for accession, suffix, repeat_sequence in repeat_records:
        if len(repeat_sequence) < config.min_repeat_length:
            continue

        signal = signals_by_accession.get(accession)
        if signal is None:
            continue

        projection = choose_projection(signal, use_trimmed=False, signal_source=config.signal_source)
        if projection is None or projection.size == 0:
            continue

        repeat_window, match_method = locate_repeat_window(
            signal.full_sequence,
            repeat_sequence,
            exact_only=config.exact_only,
        )

        if repeat_window is None:
            start_index = None
            end_index = None
            window_length = None
            score_value = None
        else:
            start_index, end_index = repeat_window
            window_length = end_index - start_index
            score_value = window_max(projection, start_index, end_index)

        scored_repeats.append(SuppliedRepeatScore(
            accession=accession,
            suffix=suffix,
            repeat_sequence=repeat_sequence,
            match_method=match_method,
            window_start=start_index,
            window_end=end_index,
            window_length=window_length,
            score=score_value,
        ))

    return scored_repeats



def write_per_accession_repeat_window_scores(repeat_scores: List[SuppliedRepeatScore], config: AnalysisConfig) -> None:
    scores_by_accession: Dict[str, List[SuppliedRepeatScore]] = {}
    for repeat_score in repeat_scores:
        scores_by_accession.setdefault(repeat_score.accession, []).append(repeat_score)

    for accession, accession_scores in scores_by_accession.items():
        output_directory = supplied_repeat_scoring_dir(config.output_dir, accession)
        output_directory.mkdir(parents=True, exist_ok=True)
        output_path = output_directory / f"{accession}_repeat_window_scores.csv"

        with output_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow([
                "accession",
                "signal_source",
                "repeat_suffix",
                "repeat_label",
                "repeat_sequence",
                "match_method",
                "window_start",
                "window_end",
                "window_length",
                "window_max_score",
            ])
            for score_record in accession_scores:
                repeat_label = f"R{score_record.suffix}" if score_record.suffix != "F" else "RF"
                writer.writerow([
                    score_record.accession,
                    config.signal_source,
                    score_record.suffix,
                    repeat_label,
                    score_record.repeat_sequence,
                    score_record.match_method,
                    "" if score_record.window_start is None else score_record.window_start,
                    "" if score_record.window_end is None else score_record.window_end,
                    "" if score_record.window_length is None else score_record.window_length,
                    "" if score_record.score is None or pd.isna(score_record.score) else round(float(score_record.score), 6),
                ])



def build_supplied_repeat_wide_table(repeat_scores: List[SuppliedRepeatScore]) -> pd.DataFrame:
    values_by_accession: Dict[str, Dict[str, float | None]] = {}

    for repeat_score in repeat_scores:
        existing_value = values_by_accession.setdefault(repeat_score.accession, {}).get(repeat_score.suffix)
        new_value = repeat_score.score
        if existing_value is None:
            values_by_accession[repeat_score.accession][repeat_score.suffix] = new_value
        elif new_value is not None and (pd.isna(existing_value) or new_value > existing_value):
            values_by_accession[repeat_score.accession][repeat_score.suffix] = new_value

    numeric_suffixes = set()
    for suffix_map in values_by_accession.values():
        for suffix in suffix_map:
            if suffix.isdigit():
                numeric_suffixes.add(int(suffix))

    maximum_repeat_number = max(numeric_suffixes) if numeric_suffixes else 0
    column_names = ["accession"] + [f"R{index}" for index in range(1, maximum_repeat_number + 1)] + ["RF"]

    rows = []
    for accession in sorted(values_by_accession.keys()):
        suffix_map = values_by_accession[accession]
        row = {"accession": accession}

        for index in range(1, maximum_repeat_number + 1):
            row[f"R{index}"] = suffix_map.get(str(index), "")

        if "F" in suffix_map and pd.notna(suffix_map["F"]):
            row["RF"] = suffix_map["F"]
        else:
            numeric_values = sorted(int(suffix) for suffix in suffix_map if suffix.isdigit())
            if numeric_values:
                highest_numeric_repeat = numeric_values[-1]
                row["RF"] = suffix_map.get(str(highest_numeric_repeat), "")
            else:
                row["RF"] = ""

        rows.append(row)

    wide_table = pd.DataFrame(rows, columns=column_names)
    for column_name in wide_table.columns:
        if column_name.startswith("R") or column_name == "RF":
            wide_table[column_name] = wide_table[column_name].apply(
                lambda value: "" if (value == "" or pd.isna(value)) else round(float(value), 6)
            )
    return wide_table



def write_supplied_repeat_summary(repeat_scores: List[SuppliedRepeatScore], config: AnalysisConfig) -> None:
    if not repeat_scores:
        return

    output_directory = global_supplied_repeat_scoring_dir(config.output_dir)
    output_directory.mkdir(parents=True, exist_ok=True)

    summary_table = build_supplied_repeat_wide_table(repeat_scores)
    summary_table.to_csv(output_directory / "supplied_repeat_scores_summary.csv", index=False)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main() -> None:
    config = parse_args()
    validate_config(config)

    config.output_dir.mkdir(parents=True, exist_ok=True)
    write_analysis_metadata(config)

    sequences_by_accession = load_full_sequences(config.input_fasta)
    trimming_map = load_trimmings(config.trimmings_csv) if config.trimmings_csv else None
    supplied_repeat_records = load_supplied_repeats(config.repeats_fasta) if config.repeats_fasta else []

    signals_by_accession: Dict[str, SequenceSignal] = {}

    for accession, full_sequence in sequences_by_accession.items():
        print(f"Processing accession {accession}")
        signal = build_full_sequence_signal(accession, full_sequence, config)

        if trimming_map is not None:
            signal = add_trimmed_signal(signal, trimming_map, config)

        signals_by_accession[accession] = signal

        if config.full_heatmaps:
            write_full_heatmap_outputs(signal, config)

        if config.scrambled_heatmaps:
            write_scrambled_heatmap_outputs(signal, config)

    if config.repeat_discovery:
        repeat_discovery_results: List[RepeatDiscoveryResult] = []
        for signal in signals_by_accession.values():
            result = run_repeat_discovery(signal, config)
            if result is None:
                continue
            repeat_discovery_results.append(result)
            write_repeat_discovery_outputs(signal, result, config)

        write_trimmed_sequences_fasta(signals_by_accession, config.output_dir)
        write_discovery_results_summary(repeat_discovery_results, signals_by_accession, config)

    if config.score_supplied_repeats:
        repeat_scores = score_supplied_repeats(signals_by_accession, supplied_repeat_records, config)
        write_per_accession_repeat_window_scores(repeat_scores, config)
        write_supplied_repeat_summary(repeat_scores, config)

    print("Done.")


if __name__ == "__main__":
    main()
