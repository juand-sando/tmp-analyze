# tapeseq

`tapeseq` is a small command-line toolkit for sequence self-similarity analysis.

All downstream analysis uses a standardized full-sequence projection computed from a reproducible shuffled baseline derived from the input sequence, accession, and master seed.

Full-sequence analysis artifacts are cached outside run outputs by default in `.tapeseq_cache/` so later modes can reuse them.

Each run also writes metadata under `output_dir/run_metadata/`, including:

- `analysis_parameters.json`
- `run.log`
- `current_state.json`

It currently exposes three modes:

- `exploratory`: write raw and smoothed heatmaps plus projection plots for each accession
- `periodicity`: run Fourier analysis on trimmed sequences and write a periodicity summary
- `repeat-analysis`: score supplied repeats against the full-sequence projection

## Installation

```bash
pip install -e .
```

## Conda Environment

From the repository root:

```bash
conda env create -f environment.yml
conda activate tapeseq
```

This installs the required packages and installs this repository in editable mode, so you can run `tapeseq` directly inside the conda environment.

## Commands

```bash
tapeseq exploratory INPUT.fasta OUTPUT_DIR --trimmings TRIMS.csv
tapeseq periodicity INPUT.fasta OUTPUT_DIR --trimmings TRIMS.csv
tapeseq repeat-analysis INPUT.fasta OUTPUT_DIR --repeats-fasta REPEATS.fasta

# Optional cache controls
tapeseq exploratory INPUT.fasta OUTPUT_DIR --cache-dir /path/to/cache
tapeseq exploratory INPUT.fasta OUTPUT_DIR --rebuild-cache
```

Verify the environment with:

```bash
tapeseq --help
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

You can also run the repository-local compatibility entrypoint:

```bash
python3 unified_sequence_analysis.py --help
```
