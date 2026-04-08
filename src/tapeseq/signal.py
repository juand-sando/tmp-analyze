from __future__ import annotations

import random
from functools import lru_cache

import numpy as np

from .models import AnalysisConfig, SequenceSignal
from .parsing import trimming_bounds


@lru_cache(maxsize=1)
def load_blosum45():
    from Bio.Align import substitution_matrices

    return substitution_matrices.load("BLOSUM45")


def score_pair(residue_a: str, residue_b: str, default: int = 0) -> int:
    matrix = load_blosum45()
    residue_a = residue_a.upper()
    residue_b = residue_b.upper()
    try:
        return matrix[(residue_a, residue_b)]
    except KeyError:
        try:
            return matrix[(residue_b, residue_a)]
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
) -> tuple[float, float]:
    randomized_projection_values: list[np.ndarray] = []

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


def add_trimmed_signal(
    signal: SequenceSignal,
    trim_map: dict[int, tuple[int, int]],
    config: AnalysisConfig,
) -> SequenceSignal:
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
