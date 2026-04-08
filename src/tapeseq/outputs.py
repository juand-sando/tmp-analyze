from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from .models import AnalysisConfig, PeriodicityResult, SequenceSignal, SuppliedRepeatScore


def per_accession_base_dir(output_dir: Path, accession: str) -> Path:
    return output_dir / "per_accession" / accession


def exploratory_dir(output_dir: Path, accession: str) -> Path:
    return per_accession_base_dir(output_dir, accession) / "exploratory"


def periodicity_dir(output_dir: Path, accession: str) -> Path:
    return per_accession_base_dir(output_dir, accession) / "periodicity"


def repeat_analysis_dir(output_dir: Path, accession: str) -> Path:
    return per_accession_base_dir(output_dir, accession) / "repeat_analysis"


def global_outputs_dir(output_dir: Path) -> Path:
    return output_dir / "global_outputs"


def global_periodicity_dir(output_dir: Path) -> Path:
    return global_outputs_dir(output_dir) / "periodicity"


def global_repeat_analysis_dir(output_dir: Path) -> Path:
    return global_outputs_dir(output_dir) / "repeat_analysis"


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


def write_wrapped_fasta_record(handle, header: str, sequence: str, line_width: int = 60) -> None:
    handle.write(f">{header}\n")
    for start_index in range(0, len(sequence), line_width):
        handle.write(sequence[start_index:start_index + line_width] + "\n")


def write_trimmed_sequences_fasta(signals_by_accession: dict[str, SequenceSignal], output_dir: Path) -> int:
    trimmed_signals = [signal for signal in signals_by_accession.values() if signal.trimmed_sequence]
    if not trimmed_signals:
        return 0

    output_directory = global_outputs_dir(output_dir)
    output_directory.mkdir(parents=True, exist_ok=True)

    with (output_directory / "trimmed_sequences.fasta").open("w", encoding="utf-8") as handle:
        for signal in trimmed_signals:
            write_wrapped_fasta_record(handle, signal.accession, signal.trimmed_sequence or "")

    return len(trimmed_signals)


def write_matrix_csv(matrix: np.ndarray, sequence: str, out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([""] + list(sequence))
        for residue, row in zip(sequence, matrix):
            writer.writerow([residue] + [f"{float(value):.4f}" for value in row])


def write_periodicity_summary(
    results: list[PeriodicityResult],
    signals_by_accession: dict[str, SequenceSignal],
    config: AnalysisConfig,
) -> None:
    if not results:
        return

    output_directory = global_periodicity_dir(config.output_dir)
    output_directory.mkdir(parents=True, exist_ok=True)

    summary_path = output_directory / "results_summary.csv"
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


def write_per_accession_repeat_window_scores(repeat_scores: list[SuppliedRepeatScore], config: AnalysisConfig) -> None:
    import pandas as pd

    scores_by_accession: dict[str, list[SuppliedRepeatScore]] = {}
    for repeat_score in repeat_scores:
        scores_by_accession.setdefault(repeat_score.accession, []).append(repeat_score)

    for accession, accession_scores in scores_by_accession.items():
        output_directory = repeat_analysis_dir(config.output_dir, accession)
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


def build_supplied_repeat_wide_table(repeat_scores: list[SuppliedRepeatScore]):
    import pandas as pd

    values_by_accession: dict[str, dict[str, float | None]] = {}

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


def write_supplied_repeat_summary(repeat_scores: list[SuppliedRepeatScore], config: AnalysisConfig) -> None:
    if not repeat_scores:
        return

    output_directory = global_repeat_analysis_dir(config.output_dir)
    output_directory.mkdir(parents=True, exist_ok=True)

    summary_table = build_supplied_repeat_wide_table(repeat_scores)
    summary_table.to_csv(output_directory / "supplied_repeat_scores_summary.csv", index=False)
