from __future__ import annotations

import csv
import re
import sys
from pathlib import Path


REPEAT_HEADER_PATTERN = re.compile(r"^(?P<accession>.+?)_R(?P<suffix>\d+|F)$")


def extract_accession(record) -> str:
    return str(record.description or record.id).split()[0]


def load_full_sequences(fasta_path: Path) -> dict[str, str]:
    from Bio import SeqIO

    sequences_by_accession: dict[str, str] = {}
    duplicate_accessions: list[str] = []
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        accession = extract_accession(record)
        if accession in sequences_by_accession:
            duplicate_accessions.append(accession)
            continue
        sequences_by_accession[accession] = str(record.seq).upper()

    if duplicate_accessions:
        duplicates = ", ".join(sorted(set(duplicate_accessions)))
        raise SystemExit(f"ERROR: duplicate accessions found in FASTA: {duplicates}")

    if not sequences_by_accession:
        raise SystemExit(f"ERROR: no sequences found in {fasta_path}")

    return sequences_by_accession


def parse_repeat_header(header: str) -> tuple[str | None, str | None]:
    token = str(header).split()[0]
    match = REPEAT_HEADER_PATTERN.match(token)
    if not match:
        return None, None
    return match.group("accession"), match.group("suffix")


def load_supplied_repeats(repeats_fasta: Path) -> list[tuple[str, str, str]]:
    from Bio import SeqIO

    repeat_records: list[tuple[str, str, str]] = []
    skipped_headers = 0
    for record in SeqIO.parse(str(repeats_fasta), "fasta"):
        accession, suffix = parse_repeat_header(str(record.description or record.id))
        if accession is None or suffix is None:
            skipped_headers += 1
            continue
        repeat_sequence = re.sub(r"[\s\-]", "", str(record.seq).upper())
        repeat_records.append((accession, suffix, repeat_sequence))

    if skipped_headers:
        print(
            f"WARNING: skipped {skipped_headers} repeat records with headers not matching ACCESSION_R<number|F>.",
            file=sys.stderr,
        )

    if not repeat_records:
        raise SystemExit(f"ERROR: no valid repeat records found in {repeats_fasta}")

    return repeat_records


def load_trimmings(csv_path: Path) -> dict[int, tuple[int, int]]:
    trimming_map: dict[int, tuple[int, int]] = {}
    skipped_rows = 0

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
                skipped_rows += 1
                continue

            if start_value < 1:
                skipped_rows += 1
                continue

            trimming_map[length_value] = (start_value, end_value)

    if skipped_rows:
        print(f"WARNING: skipped {skipped_rows} invalid trimming rows in {csv_path}.", file=sys.stderr)

    if not trimming_map:
        raise SystemExit(f"ERROR: no valid trimming rows found in {csv_path}.")

    return trimming_map


def trimming_bounds(full_length: int, trim_map: dict[int, tuple[int, int]]) -> tuple[int, int] | None:
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
