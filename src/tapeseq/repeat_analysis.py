from __future__ import annotations

import numpy as np

from .models import AnalysisConfig, SequenceSignal, SuppliedRepeatScore
from .outputs import (
    write_analysis_metadata,
    write_per_accession_repeat_window_scores,
    write_supplied_repeat_summary,
)
from .parsing import load_full_sequences, load_supplied_repeats, load_trimmings
from .signal import add_trimmed_signal, build_or_load_full_sequence_signal, choose_projection


def find_repeat_range_exact(full_sequence: str, repeat_sequence: str) -> tuple[int, int] | None:
    start_index = full_sequence.find(repeat_sequence)
    if start_index < 0:
        return None
    return start_index, start_index + len(repeat_sequence)


def find_repeat_range_align(full_sequence: str, repeat_sequence: str) -> tuple[int, int] | None:
    from Bio import pairwise2

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


def locate_repeat_window(
    full_sequence: str,
    repeat_sequence: str,
    exact_only: bool = False,
) -> tuple[tuple[int, int] | None, str]:
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
    signals_by_accession: dict[str, SequenceSignal],
    repeat_records: list[tuple[str, str, str]],
    config: AnalysisConfig,
) -> list[SuppliedRepeatScore]:
    scored_repeats: list[SuppliedRepeatScore] = []

    for accession, suffix, repeat_sequence in repeat_records:
        if len(repeat_sequence) < config.min_repeat_length:
            continue

        signal = signals_by_accession.get(accession)
        if signal is None:
            continue

        projection = choose_projection(signal, use_trimmed=False)
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

        scored_repeats.append(
            SuppliedRepeatScore(
                accession=accession,
                suffix=suffix,
                repeat_sequence=repeat_sequence,
                match_method=match_method,
                window_start=start_index,
                window_end=end_index,
                window_length=window_length,
                score=score_value,
            )
        )

    return scored_repeats


def run(config: AnalysisConfig) -> None:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    write_analysis_metadata(config)

    sequences_by_accession = load_full_sequences(config.input_fasta)
    trimming_map = load_trimmings(config.trimmings_csv) if config.trimmings_csv else None
    repeat_records = load_supplied_repeats(config.repeats_fasta) if config.repeats_fasta else []

    signals_by_accession: dict[str, SequenceSignal] = {}
    for accession, full_sequence in sequences_by_accession.items():
        print(f"Processing accession {accession}")
        signal = build_or_load_full_sequence_signal(accession, full_sequence, config)
        if trimming_map is not None:
            signal = add_trimmed_signal(signal, trimming_map, config)
        signals_by_accession[accession] = signal

    repeat_scores = score_supplied_repeats(signals_by_accession, repeat_records, config)
    write_per_accession_repeat_window_scores(repeat_scores, config)
    write_supplied_repeat_summary(repeat_scores, config)

    print(
        f"Done. Processed {len(signals_by_accession)} accession(s); "
        f"scored {len(repeat_scores)} supplied repeat record(s)."
    )
