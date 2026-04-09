# tapemeasure-analyze

This repository now installs two command-line tools for tape measure protein analysis:

- `tapeseq`: sequence-focused self-similarity, periodicity, repeat, and domain-mapping analysis
- `tapebulk`: structure-coordinate bulkiness and pairwise-distance analysis

Both tools live in the same Python distribution and share the same installation.

## `tapeseq`

`tapeseq` is the sequence-analysis toolkit.

All downstream analysis uses a standardized full-sequence projection computed from a reproducible shuffled baseline derived from the input sequence, accession, and master seed.

Full-sequence analysis artifacts are cached outside run outputs by default in `.tapeseq_cache/` so later modes can reuse them.

Each run also writes metadata under `output_dir/run_metadata/`, including:

- `analysis_parameters.json`
- `run.log`
- `current_state.json`

It exposes four modes:

- `exploratory`: write raw and smoothed heatmaps plus projection plots for each accession
- `periodicity`: run Fourier analysis on the full or trimmed focused window and write a periodicity summary
- `repeat-analysis`: score supplied repeats against the full-sequence projection
- `map-domain`: project a reference domain through an alignment and generate accession-based trimming rules

## Installation

```bash
pip install -e .
```

After updating an existing local install from this repo, run `pip install -e .` again so the new `tapebulk` console script is installed.

## Conda Environment

From the repository root:

```bash
conda env create -f environment.yml
conda activate tapeseq
```

This installs the required packages and installs this repository in editable mode, so you can run both `tapeseq` and `tapebulk` inside the conda environment.

## `tapeseq` Commands

```bash
tapeseq exploratory -i INPUT.fasta -o OUTPUT_DIR --trimmings TRIMS.csv
tapeseq periodicity -i INPUT.fasta -o OUTPUT_DIR --trimmings TRIMS.csv
tapeseq repeat-analysis -i INPUT.fasta -o OUTPUT_DIR --repeats-fasta REPEATS.fasta
tapeseq map-domain -a ALIGNMENT.fasta -o OUTPUT_DIR --reference-accession REF_ACC --reference-domain DOMAIN.fasta

# Optional cache controls
tapeseq exploratory -i INPUT.fasta -o OUTPUT_DIR --cache-dir /path/to/cache
tapeseq exploratory -i INPUT.fasta -o OUTPUT_DIR --rebuild-cache
```

Verify the environment with:

```bash
tapeseq --help
tapebulk --help
```

## Trimming CSV

The trimming CSV must contain `start` and `end`, plus either `tmp_length` or `accession`.

Length-based rules:

```csv
tmp_length,start,end
120,10,95
150,20,130
```

Accession-based rules:

```csv
accession,start,end
ACC001,10,95
ACC002,20,130
```

Mixed files are also accepted. When both are present, accession-specific rules take precedence over length-based rules.

## Domain Mapping

`map-domain` takes an aligned FASTA, a reference accession from that alignment, and a FASTA containing the reference domain sequence. It projects that domain window through the alignment and writes:

- `global_outputs/map_domain/mapped_domain_trimmings.csv`
- `global_outputs/map_domain/mapped_domain_qc.csv`

Only conservative accession-based trimming rules are emitted. Accessions with weak or ambiguous mappings are kept in the QC file and skipped from the trimming CSV.

You can also run the repository-local compatibility entrypoint:

```bash
python3 unified_sequence_analysis.py --help
```

## `tapebulk`

`tapebulk` is the structure-coordinate toolkit. It currently takes a CSV with columns `Pos`, `Res`, `X`, `Y`, `Z`.

It exposes two modes:

- `bulkiness`: neighbor counts, optional straight-chain baseline subtraction, and kNN-derived bulkiness scores
- `pairwise-histogram`: all pairwise distances plus a histogram

Example commands:

```bash
tapebulk bulkiness -i POINTS.csv -o OUTPUT_DIR --radii 6 10 15 --k-neighbors 10 --subtract-baseline --make-plots
tapebulk pairwise-histogram -i POINTS.csv -o OUTPUT_DIR --histogram-bins 40
```

`tapebulk` writes its outputs under `OUTPUT_DIR/global_outputs/bulkiness/` plus the shared `run_metadata/` directory.

The repository-local compatibility wrapper for the older standalone script is still available:

```bash
python3 bulkiness_v2.py --help
```
