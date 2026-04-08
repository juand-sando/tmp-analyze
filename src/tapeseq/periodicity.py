from __future__ import annotations

import math

import numpy as np

from .models import AnalysisConfig, PeriodicityResult, SequenceSignal
from .outputs import (
    periodicity_dir,
    write_analysis_metadata,
    write_matrix_csv,
    write_periodicity_summary,
    write_trimmed_sequences_fasta,
    write_wrapped_fasta_record,
)
from .parsing import load_full_sequences, load_trimmings
from .plotting import save_fft_plot, save_projection_plot
from .signal import add_trimmed_signal, build_full_sequence_signal, choose_projection


def fft_power_and_peak(signal: np.ndarray, minimum_frequency: int, maximum_frequency: int) -> tuple[np.ndarray, int, float]:
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


def slice_repeats(sequence: str, period: float) -> list[str]:
    if not sequence or not np.isfinite(period) or period < 1.0:
        return []

    window_size = max(1, int(math.floor(period)))
    repeat_sequences: list[str] = []
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


def analyze_periodicity(signal: SequenceSignal, config: AnalysisConfig) -> PeriodicityResult | None:
    if signal.trimmed_sequence is None:
        return None

    projection = choose_projection(signal, use_trimmed=True)
    if projection is None or projection.size == 0:
        return None

    power_spectrum, best_frequency, best_power = fft_power_and_peak(
        projection,
        minimum_frequency=config.fft_min,
        maximum_frequency=config.fft_max,
    )
    period = period_from_frequency(len(projection), best_frequency)
    repeat_sequences = slice_repeats(signal.trimmed_sequence, period)

    return PeriodicityResult(
        accession=signal.accession,
        trimmed_length=len(signal.trimmed_sequence),
        best_frequency=best_frequency,
        best_power=best_power,
        period=period,
        repeat_count=len(repeat_sequences),
        repeat_sequences=repeat_sequences,
    )


def write_periodicity_outputs(signal: SequenceSignal, result: PeriodicityResult, config: AnalysisConfig) -> None:
    output_directory = periodicity_dir(config.output_dir, signal.accession)
    output_directory.mkdir(parents=True, exist_ok=True)

    if signal.trimmed_sequence is not None:
        with (output_directory / f"{signal.accession}_trimmed.fasta").open("w", encoding="utf-8") as handle:
            write_wrapped_fasta_record(handle, signal.accession, signal.trimmed_sequence)

    if signal.smoothed_matrix_trimmed is not None and signal.trimmed_sequence is not None:
        write_matrix_csv(
            signal.smoothed_matrix_trimmed,
            signal.trimmed_sequence,
            output_directory / f"{signal.accession}_trimmed_matrix.csv",
        )

    projection = choose_projection(signal, use_trimmed=True)
    if projection is not None:
        save_projection_plot(
            projection,
            f"{signal.accession} trimmed standardized projection",
            output_directory / f"{signal.accession}_projection.png",
        )
        power_spectrum, _, _ = fft_power_and_peak(
            projection,
            minimum_frequency=config.fft_min,
            maximum_frequency=config.fft_max,
        )
        save_fft_plot(
            power_spectrum,
            f"{signal.accession} periodicity FFT power",
            output_directory / f"{signal.accession}_fft.png",
            minimum_frequency=config.fft_min,
            maximum_frequency=config.fft_max,
            highlight_frequency=result.best_frequency,
        )

    with (output_directory / f"{signal.accession}_periodicity_segments.fasta").open("w", encoding="utf-8") as handle:
        for repeat_index, repeat_sequence in enumerate(result.repeat_sequences, start=1):
            suffix = f"R{repeat_index}" if repeat_index < len(result.repeat_sequences) else "RF"
            write_wrapped_fasta_record(handle, f"{signal.accession}_{suffix}", repeat_sequence)


def run(config: AnalysisConfig) -> None:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    write_analysis_metadata(config)

    sequences_by_accession = load_full_sequences(config.input_fasta)
    trimming_map = load_trimmings(config.trimmings_csv) if config.trimmings_csv else None

    signals_by_accession: dict[str, SequenceSignal] = {}
    results: list[PeriodicityResult] = []

    for accession, full_sequence in sequences_by_accession.items():
        print(f"Processing accession {accession}")
        signal = build_full_sequence_signal(accession, full_sequence, config)
        if trimming_map is not None:
            signal = add_trimmed_signal(signal, trimming_map, config)
        signals_by_accession[accession] = signal

        result = analyze_periodicity(signal, config)
        if result is None:
            continue

        results.append(result)
        write_periodicity_outputs(signal, result, config)

    trimmed_written = write_trimmed_sequences_fasta(signals_by_accession, config.output_dir)
    write_periodicity_summary(results, signals_by_accession, config)

    print(
        f"Done. Processed {len(signals_by_accession)} accession(s); "
        f"trimmed {trimmed_written}; periodicity results for {len(results)} accession(s)."
    )
