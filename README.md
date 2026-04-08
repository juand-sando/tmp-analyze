# tapeseq

`tapeseq` is a small command-line toolkit for sequence self-similarity analysis.

All downstream analysis uses a standardized full-sequence projection computed from a reproducible shuffled baseline derived from the input sequence, accession, and master seed.

Full-sequence analysis artifacts are cached outside run outputs by default in `.tapeseq_cache/` so later modes can reuse them.

It currently exposes three modes:

- `exploratory`: write raw and smoothed heatmaps plus projection plots for each accession
- `periodicity`: run Fourier analysis on trimmed sequences and write a periodicity summary
- `repeat-analysis`: score supplied repeats against the full-sequence projection

## Installation

```bash
pip install -e .
```

## Commands

```bash
tapeseq exploratory INPUT.fasta OUTPUT_DIR --trimmings TRIMS.csv
tapeseq periodicity INPUT.fasta OUTPUT_DIR --trimmings TRIMS.csv
tapeseq repeat-analysis INPUT.fasta OUTPUT_DIR --repeats-fasta REPEATS.fasta

# Optional cache controls
tapeseq exploratory INPUT.fasta OUTPUT_DIR --cache-dir /path/to/cache
tapeseq exploratory INPUT.fasta OUTPUT_DIR --rebuild-cache
```

You can also run the repository-local compatibility entrypoint:

```bash
python3 unified_sequence_analysis.py --help
```
