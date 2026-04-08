from __future__ import annotations

import math
import sys

import numpy as np

from .models import AnalysisConfig, PeriodicityResult, SequenceSignal
from .outputs import (
    periodicity_dir,
    write_analysis_metadata,
    write_global_periodicity_segments_fasta,
    write_matrix_csv,
    write_periodicity_notice,
    write_periodicity_summary,
    write_trimmed_sequences_fasta,
    write_wrapped_fasta_record,
)
from .parsing import load_full_sequences, load_trimmings
from .plotting import save_fft_plot, save_projection_plot
from .signal import add_trimmed_signal, build_or_load_full_sequence_signal, choose_projection


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
    use_trimmed = config.trimmings_csv is not None
    analysis_sequence = signal.trimmed_sequence if use_trimmed else signal.full_sequence
    if analysis_sequence is None:
        return None

    projection = choose_projection(signal, use_trimmed=use_trimmed)
    if projection is None or projection.size == 0:
        return None

    power_spectrum, best_frequency, best_power = fft_power_and_peak(
        projection,
        minimum_frequency=config.fft_min,
        maximum_frequency=config.fft_max,
    )
    period = period_from_frequency(len(projection), best_frequency)
    repeat_sequences = slice_repeats(analysis_sequence, period)

    return PeriodicityResult(
        accession=signal.accession,
        trimmed_length=len(analysis_sequence),
        best_frequency=best_frequency,
        best_power=best_power,
        period=period,
        repeat_count=len(repeat_sequences),
        repeat_sequences=repeat_sequences,
    )


def write_periodicity_outputs(signal: SequenceSignal, result: PeriodicityResult, config: AnalysisConfig) -> None:
    output_directory = periodicity_dir(config.output_dir, signal.accession)
    output_directory.mkdir(parents=True, exist_ok=True)

    use_trimmed = config.trimmings_csv is not None and signal.trimmed_sequence is not None
    analysis_sequence = signal.trimmed_sequence if use_trimmed else signal.full_sequence
    analysis_matrix = signal.smoothed_matrix_trimmed if use_trimmed else signal.smoothed_matrix_full
    analysis_label = "trimmed" if use_trimmed else "full"

    if use_trimmed:
        with (output_directory / f"{signal.accession}_trimmed.fasta").open("w", encoding="utf-8") as handle:
            write_wrapped_fasta_record(handle, signal.accession, signal.trimmed_sequence)

    if analysis_matrix is not None and analysis_sequence is not None:
        write_matrix_csv(
            analysis_matrix,
            analysis_sequence,
            output_directory / f"{signal.accession}_{analysis_label}_matrix.csv",
        )

    projection = choose_projection(signal, use_trimmed=use_trimmed)
    if projection is not None:
        save_projection_plot(
            projection,
            f"{signal.accession} {analysis_label} standardized projection",
            output_directory / f"{signal.accession}_projection.png",
        )
        power_spectrum, _, _ = fft_power_and_peak(
            projection,
            minimum_frequency=config.fft_min,
            maximum_frequency=config.fft_max,
        )
        save_fft_plot(
            power_spectrum,
            f"{signal.accession} {analysis_label} periodicity FFT power",
            output_directory / f"{signal.accession}_fft.png",
            minimum_frequency=config.fft_min,
            maximum_frequency=config.fft_max,
            highlight_frequency=result.best_frequency,
        )

    with (output_directory / f"{signal.accession}_periodicity_segments.fasta").open("w", encoding="utf-8") as handle:
        for repeat_index, repeat_sequence in enumerate(result.repeat_sequences, start=1):
            suffix = f"R{repeat_index}" if repeat_index < len(result.repeat_sequences) else "RF"
            write_wrapped_fasta_record(handle, f"{signal.accession}_{suffix}", repeat_sequence)


def period_is_integer_like(period: float, tolerance: float = 1e-6) -> bool:
    if not np.isfinite(period):
        return False
    return abs(period - round(period)) <= tolerance


def emit_periodicity_warnings(results: list[PeriodicityResult]) -> int:
    if not results:
        return 0

    finite_periods = np.array([result.period for result in results if np.isfinite(result.period)], dtype=float)
    if finite_periods.size == 0:
        for result in results:
            print(
                f"WARNING: accession {result.accession} has a non-finite inferred period; putative repeats may be unreliable.",
                file=sys.stderr,
            )
        return len(results)

    mean_period = float(np.mean(finite_periods))
    std_period = float(np.std(finite_periods))
    warning_count = 0

    for result in results:
        warnings_for_result: list[str] = []
        if not np.isfinite(result.period):
            warnings_for_result.append("period length is non-finite, so the putative repeat segmentation may be unreliable")
        elif not period_is_integer_like(result.period):
            warnings_for_result.append(
                f"period length {result.period:.6f} is not an integer, so the putative repeat segmentation may be far from true boundaries"
            )
        if np.isfinite(result.period) and std_period > 0.0 and abs(result.period - mean_period) > std_period:
            warnings_for_result.append(
                f"period length {result.period:.6f} is more than 1 SD from the analyzed-set mean ({mean_period:.6f} +/- {std_period:.6f})"
            )

        if warnings_for_result:
            warning_count += len(warnings_for_result)
            for warning_text in warnings_for_result:
                print(f"WARNING: accession {result.accession}: {warning_text}.", file=sys.stderr)

    return warning_count


def run(config: AnalysisConfig) -> None:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    write_analysis_metadata(config)

    sequences_by_accession = load_full_sequences(config.input_fasta)
    trimming_map = load_trimmings(config.trimmings_csv) if config.trimmings_csv else None

    signals_by_accession: dict[str, SequenceSignal] = {}
    results: list[PeriodicityResult] = []

    for accession, full_sequence in sequences_by_accession.items():
        print(f"Processing accession {accession}")
        signal = build_or_load_full_sequence_signal(accession, full_sequence, config)
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
    global_repeat_records = write_global_periodicity_segments_fasta(results, config.output_dir)
    warning_count = emit_periodicity_warnings(results)
    notice_message = (
        "Putative repeats in this periodicity analysis are derived from FFT-estimated period lengths and sequence segmentation. "
        "They are not validated biological repeats and should be corroborated with additional analysis."
    )
    write_periodicity_notice(config.output_dir, notice_message)
    print(notice_message, file=sys.stderr)

    print(
        f"Done. Processed {len(signals_by_accession)} accession(s); "
        f"trimmed {trimmed_written}; periodicity results for {len(results)} accession(s); "
        f"wrote {global_repeat_records} putative repeat record(s); emitted {warning_count} warning(s)."
    )
