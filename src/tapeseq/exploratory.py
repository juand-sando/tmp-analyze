from __future__ import annotations

from .models import AnalysisConfig, SequenceSignal
from .outputs import (
    exploratory_dir,
    write_analysis_metadata,
    write_trimmed_sequences_fasta,
    write_wrapped_fasta_record,
)
from .parsing import load_full_sequences, load_trimmings
from .plotting import save_heatmap, save_projection_plot
from .signal import add_trimmed_signal, build_full_sequence_signal, choose_projection


def write_exploratory_outputs(signal: SequenceSignal, config: AnalysisConfig) -> None:
    output_directory = exploratory_dir(config.output_dir, signal.accession)
    output_directory.mkdir(parents=True, exist_ok=True)

    if signal.raw_matrix_full is not None:
        save_heatmap(
            signal.raw_matrix_full,
            f"{signal.accession} raw matrix (full)",
            output_directory / f"{signal.accession}_full_raw.png",
        )

    if signal.smoothed_matrix_full is not None:
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

    full_projection = choose_projection(signal, use_trimmed=False, signal_source=config.signal_source)
    if full_projection is not None:
        save_projection_plot(
            full_projection,
            f"{signal.accession} full projection",
            output_directory / f"{signal.accession}_full_projection.png",
        )

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

    trimmed_projection = choose_projection(signal, use_trimmed=True, signal_source=config.signal_source)
    if trimmed_projection is not None:
        save_projection_plot(
            trimmed_projection,
            f"{signal.accession} trimmed projection",
            output_directory / f"{signal.accession}_trimmed_projection.png",
        )


def run(config: AnalysisConfig) -> None:
    config.output_dir.mkdir(parents=True, exist_ok=True)
    write_analysis_metadata(config)

    sequences_by_accession = load_full_sequences(config.input_fasta)
    trimming_map = load_trimmings(config.trimmings_csv) if config.trimmings_csv else None

    signals_by_accession: dict[str, SequenceSignal] = {}
    for accession, full_sequence in sequences_by_accession.items():
        print(f"Processing accession {accession}")
        signal = build_full_sequence_signal(accession, full_sequence, config)
        if trimming_map is not None:
            signal = add_trimmed_signal(signal, trimming_map, config)
        signals_by_accession[accession] = signal
        write_exploratory_outputs(signal, config)

    trimmed_written = 0
    if trimming_map is not None:
        trimmed_written = write_trimmed_sequences_fasta(signals_by_accession, config.output_dir)

    print(
        f"Done. Processed {len(signals_by_accession)} accession(s); "
        f"wrote trimmed sequences for {trimmed_written} accession(s)."
    )
