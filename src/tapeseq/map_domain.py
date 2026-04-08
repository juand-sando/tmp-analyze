from __future__ import annotations

import csv
import sys
from dataclasses import asdict

from .models import AnalysisConfig, MapDomainQCRecord
from .outputs import global_outputs_dir, write_analysis_metadata
from .parsing import extract_accession
from .run_logging import get_logger, log_file_write


GAP_CHARACTERS = {"-", "."}
MIN_REFERENCE_ALIGNMENT_COVERAGE = 0.90
MIN_REFERENCE_ALIGNMENT_IDENTITY = 0.90
MIN_PROJECTED_DOMAIN_COVERAGE = 0.80


def ungap(sequence: str) -> str:
    return "".join(character for character in sequence if character not in GAP_CHARACTERS)


def load_alignment_sequences(fasta_path):
    from Bio import SeqIO

    sequences_by_accession: dict[str, str] = {}
    duplicate_accessions: list[str] = []
    alignment_lengths: set[int] = set()

    for record in SeqIO.parse(str(fasta_path), "fasta"):
        accession = extract_accession(record)
        aligned_sequence = str(record.seq).upper()
        if accession in sequences_by_accession:
            duplicate_accessions.append(accession)
            continue
        sequences_by_accession[accession] = aligned_sequence
        alignment_lengths.add(len(aligned_sequence))

    if duplicate_accessions:
        duplicates = ", ".join(sorted(set(duplicate_accessions)))
        raise SystemExit(f"ERROR: duplicate accessions found in alignment FASTA: {duplicates}")

    if not sequences_by_accession:
        raise SystemExit(f"ERROR: no sequences found in alignment FASTA {fasta_path}")

    if len(alignment_lengths) != 1:
        raise SystemExit("ERROR: alignment FASTA sequences do not all have the same aligned length.")

    return sequences_by_accession


def load_single_fasta_sequence(fasta_path):
    from Bio import SeqIO

    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(records) != 1:
        raise SystemExit(f"ERROR: expected exactly one sequence in {fasta_path}, found {len(records)}.")
    record = records[0]
    return extract_accession(record), str(record.seq).upper()


def locate_domain_exact(reference_sequence: str, domain_sequence: str):
    start_index = reference_sequence.find(domain_sequence)
    if start_index < 0:
        return None
    end_index = start_index + len(domain_sequence)
    return {
        "method": "exact",
        "start": start_index,
        "end": end_index,
        "coverage": 1.0,
        "identity": 1.0,
        "warning": "",
    }


def locate_domain_alignment(reference_sequence: str, domain_sequence: str):
    from Bio import pairwise2

    alignments = pairwise2.align.localms(reference_sequence, domain_sequence, 2, -1, -5, -1)
    if not alignments:
        return None

    best_alignment = alignments[0]
    if hasattr(best_alignment, "seqA"):
        aligned_reference = best_alignment.seqA
        aligned_domain = best_alignment.seqB
        reference_start = int(best_alignment.start)
    else:
        aligned_reference, aligned_domain, _score, reference_start, _end = best_alignment

    reference_position = reference_start
    mapped_start = None
    mapped_end = None
    domain_residue_count = len(domain_sequence)
    aligned_domain_residues = 0
    identical_residues = 0

    for reference_character, domain_character in zip(aligned_reference, aligned_domain):
        current_reference_position = None
        if reference_character != "-":
            current_reference_position = reference_position
            reference_position += 1

        if domain_character == "-":
            continue

        if current_reference_position is not None:
            if mapped_start is None:
                mapped_start = current_reference_position
            mapped_end = current_reference_position
            aligned_domain_residues += 1
            if reference_character == domain_character:
                identical_residues += 1

    if mapped_start is None or mapped_end is None or domain_residue_count == 0:
        return None

    coverage = aligned_domain_residues / domain_residue_count
    identity = identical_residues / aligned_domain_residues if aligned_domain_residues else 0.0
    return {
        "method": "alignment",
        "start": mapped_start,
        "end": mapped_end + 1,
        "coverage": coverage,
        "identity": identity,
        "warning": (
            f"reference domain matched by local alignment with coverage {coverage:.3f} "
            f"and identity {identity:.3f}"
        ),
    }


def locate_domain_in_reference(reference_sequence: str, domain_sequence: str):
    exact = locate_domain_exact(reference_sequence, domain_sequence)
    if exact is not None:
        return exact

    aligned = locate_domain_alignment(reference_sequence, domain_sequence)
    if aligned is None:
        raise SystemExit("ERROR: reference domain sequence could not be matched to the reference accession sequence.")

    if (
        aligned["coverage"] < MIN_REFERENCE_ALIGNMENT_COVERAGE
        or aligned["identity"] < MIN_REFERENCE_ALIGNMENT_IDENTITY
    ):
        raise SystemExit(
            "ERROR: reference domain alignment to the reference accession was too weak "
            f"(coverage={aligned['coverage']:.3f}, identity={aligned['identity']:.3f})."
        )

    print(f"WARNING: {aligned['warning']}.", file=sys.stderr)
    get_logger().warning(f"warning: {aligned['warning']}.")
    return aligned


def residue_index_to_alignment_column_map(aligned_sequence: str) -> list[int]:
    return [column_index for column_index, character in enumerate(aligned_sequence) if character not in GAP_CHARACTERS]


def reference_window_columns(aligned_reference_sequence: str, domain_start: int, domain_end: int) -> tuple[int, int]:
    residue_columns = residue_index_to_alignment_column_map(aligned_reference_sequence)
    if domain_end <= domain_start or domain_end > len(residue_columns):
        raise SystemExit("ERROR: mapped reference domain bounds are out of range for the aligned reference sequence.")
    start_column = residue_columns[domain_start]
    end_column = residue_columns[domain_end - 1] + 1
    return start_column, end_column


def map_alignment_window_to_accession(
    accession: str,
    aligned_sequence: str,
    aligned_reference_window: str,
    start_column: int,
    end_column: int,
    reference_accession: str,
    reference_domain_method: str,
    reference_domain_start: int,
    reference_domain_end: int,
    reference_domain_length: int,
) -> MapDomainQCRecord:
    aligned_window = aligned_sequence[start_column:end_column]
    non_gap_columns = [index for index, character in enumerate(aligned_window) if character not in GAP_CHARACTERS]

    if not non_gap_columns:
        return MapDomainQCRecord(
            accession=accession,
            status="unmapped",
            warning="projected alignment window contains only gaps for this accession",
            reference_accession=reference_accession,
            reference_domain_method=reference_domain_method,
            reference_domain_start=reference_domain_start,
            reference_domain_end=reference_domain_end,
            alignment_start_column=start_column,
            alignment_end_column=end_column,
            reference_domain_length=reference_domain_length,
            mapped_start=None,
            mapped_end=None,
            mapped_length=None,
            coverage_fraction=0.0,
            edge_gap_left=len(aligned_window),
            edge_gap_right=len(aligned_window),
        )

    first_non_gap_offset = non_gap_columns[0]
    last_non_gap_offset = non_gap_columns[-1]
    first_non_gap_column = start_column + first_non_gap_offset
    last_non_gap_column = start_column + last_non_gap_offset

    mapped_start_zero_based = sum(
        1 for character in aligned_sequence[:first_non_gap_column] if character not in GAP_CHARACTERS
    )
    mapped_end_exclusive = sum(
        1 for character in aligned_sequence[:last_non_gap_column + 1] if character not in GAP_CHARACTERS
    )
    mapped_length = mapped_end_exclusive - mapped_start_zero_based

    reference_domain_columns = [
        index
        for index, character in enumerate(aligned_reference_window)
        if character not in GAP_CHARACTERS
    ]
    supported_reference_positions = sum(
        1
        for column_offset in reference_domain_columns
        if aligned_window[column_offset] not in GAP_CHARACTERS
    )
    coverage_fraction = (
        supported_reference_positions / reference_domain_length if reference_domain_length > 0 else 0.0
    )

    edge_gap_left = first_non_gap_offset
    edge_gap_right = len(aligned_window) - last_non_gap_offset - 1
    warnings: list[str] = []

    if coverage_fraction < MIN_PROJECTED_DOMAIN_COVERAGE:
        warnings.append(
            f"coverage {coverage_fraction:.3f} is below the conservative threshold {MIN_PROJECTED_DOMAIN_COVERAGE:.2f}"
        )
    if edge_gap_left > 0:
        warnings.append(f"left edge begins with {edge_gap_left} gap column(s)")
    if edge_gap_right > 0:
        warnings.append(f"right edge ends with {edge_gap_right} gap column(s)")

    status = "ok" if not warnings else "partial"
    mapped_start = mapped_start_zero_based + 1
    mapped_end = mapped_end_exclusive

    return MapDomainQCRecord(
        accession=accession,
        status=status,
        warning="; ".join(warnings),
        reference_accession=reference_accession,
        reference_domain_method=reference_domain_method,
        reference_domain_start=reference_domain_start,
        reference_domain_end=reference_domain_end,
        alignment_start_column=start_column,
        alignment_end_column=end_column,
        reference_domain_length=reference_domain_length,
        mapped_start=mapped_start,
        mapped_end=mapped_end,
        mapped_length=mapped_length,
        coverage_fraction=coverage_fraction,
        edge_gap_left=edge_gap_left,
        edge_gap_right=edge_gap_right,
    )


def write_map_domain_outputs(config: AnalysisConfig, qc_records: list[MapDomainQCRecord]) -> tuple[int, int]:
    output_directory = global_outputs_dir(config.output_dir) / "map_domain"
    output_directory.mkdir(parents=True, exist_ok=True)

    trimmings_path = output_directory / "mapped_domain_trimmings.csv"
    qc_path = output_directory / "mapped_domain_qc.csv"

    trimming_records = [record for record in qc_records if record.status == "ok" and record.mapped_start is not None and record.mapped_end is not None]

    log_file_write(trimmings_path)
    with trimmings_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["accession", "start", "end", "reference_accession", "reference_domain_method"])
        for record in trimming_records:
            writer.writerow([
                record.accession,
                record.mapped_start,
                record.mapped_end,
                record.reference_accession,
                record.reference_domain_method,
            ])

    log_file_write(qc_path)
    with qc_path.open("w", newline="", encoding="utf-8") as handle:
        fieldnames = list(asdict(qc_records[0]).keys()) if qc_records else list(asdict(
            MapDomainQCRecord(
                accession="",
                status="",
                warning="",
                reference_accession="",
                reference_domain_method="",
                reference_domain_start=None,
                reference_domain_end=None,
                alignment_start_column=None,
                alignment_end_column=None,
                reference_domain_length=0,
                mapped_start=None,
                mapped_end=None,
                mapped_length=None,
                coverage_fraction=None,
                edge_gap_left=None,
                edge_gap_right=None,
            )
        ).keys())
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in qc_records:
            writer.writerow(asdict(record))

    return len(trimming_records), len(qc_records)


def run(config: AnalysisConfig) -> str:
    write_analysis_metadata(config)

    alignment_by_accession = load_alignment_sequences(config.input_fasta)
    if config.reference_accession not in alignment_by_accession:
        raise SystemExit(
            f"ERROR: reference accession {config.reference_accession} is not present in the alignment FASTA."
        )

    _domain_header, reference_domain_sequence = load_single_fasta_sequence(config.reference_domain_fasta)
    reference_domain_sequence = ungap(reference_domain_sequence.upper())
    if not reference_domain_sequence:
        raise SystemExit("ERROR: reference domain FASTA does not contain any residues after removing gaps.")

    aligned_reference_sequence = alignment_by_accession[config.reference_accession]
    reference_sequence = ungap(aligned_reference_sequence)
    reference_match = locate_domain_in_reference(reference_sequence, reference_domain_sequence)
    start_column, end_column = reference_window_columns(
        aligned_reference_sequence,
        reference_match["start"],
        reference_match["end"],
    )
    aligned_reference_window = aligned_reference_sequence[start_column:end_column]

    qc_records: list[MapDomainQCRecord] = []
    skipped_records = 0
    for accession, aligned_sequence in alignment_by_accession.items():
        qc_record = map_alignment_window_to_accession(
            accession=accession,
            aligned_sequence=aligned_sequence,
            aligned_reference_window=aligned_reference_window,
            start_column=start_column,
            end_column=end_column,
            reference_accession=config.reference_accession,
            reference_domain_method=reference_match["method"],
            reference_domain_start=reference_match["start"] + 1,
            reference_domain_end=reference_match["end"],
            reference_domain_length=len(reference_domain_sequence),
        )
        qc_records.append(qc_record)
        if qc_record.status != "ok":
            skipped_records += 1
            warning_suffix = f": {qc_record.warning}" if qc_record.warning else ""
            print(
                f"WARNING: accession {accession} was not emitted as a trimming rule ({qc_record.status}){warning_suffix}.",
                file=sys.stderr,
            )
            get_logger().warning(
                f"warning: accession {accession} was not emitted as a trimming rule ({qc_record.status}){warning_suffix}."
            )

    trimming_count, qc_count = write_map_domain_outputs(config, qc_records)
    summary = (
        f"Done. Processed {len(alignment_by_accession)} aligned accession(s); "
        f"wrote {trimming_count} conservative trimming rule(s); "
        f"recorded QC for {qc_count} accession(s); "
        f"skipped {skipped_records} accession(s)."
    )
    print(summary)
    return summary
