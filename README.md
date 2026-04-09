# tapemeasure-analyze

`tapemeasure-analyze` is a small research toolkit for studying tape measure proteins with two complementary command-line tools:

- `tapeseq` for sequence-based analysis
- `tapebulk` for structure-coordinate analysis

Both tools are installed together and are meant to live in the same conda environment.

**Biology**
Tape measure proteins are typically very long phage proteins whose lengths and internal organization are often linked to tail architecture and assembly. In practice, people studying them often want to ask several related questions:

- Does the sequence contain repeated or quasi-repeated internal structure?
- Are there periodic patterns that may reflect recurring modules?
- Can a domain annotated in one protein be mapped conservatively onto homologous proteins?
- Does a structural model show locally crowded or expanded regions along the chain?

This repository is built around those questions. The tools here are intended to help generate and compare hypotheses across tape measure proteins; they do not by themselves prove structural or functional mechanisms.

**What This Repository Contains**

- `tapeseq`: computes sequence self-similarity signals, periodicity summaries, repeat-associated scores, and alignment-based domain mappings
- `tapebulk`: computes local coordinate crowding and pairwise distance distributions from 3D residue coordinates

Internally, the code is split into:

- [`src/tapemeasure/seq`](/Users/juand.sando/Documents/fv3_repositories/tapeseq-comp/src/tapemeasure/seq) for sequence workflows
- [`src/tapemeasure/bulk`](/Users/juand.sando/Documents/fv3_repositories/tapeseq-comp/src/tapemeasure/bulk) for structure workflows
- [`src/tapemeasure/shared`](/Users/juand.sando/Documents/fv3_repositories/tapeseq-comp/src/tapemeasure/shared) for shared run metadata/logging helpers

## Installation

**Who This Is For**

These instructions assume you are new to Python environments and just want a working install on your computer.

**Prerequisite**

Install Miniconda or Anaconda if you do not already have conda available in your terminal.

Check that conda exists:

```bash
conda --version
```

If that prints a version number, continue.

**Clone The Repository**

```bash
git clone <YOUR-REPO-URL>
cd tapeseq-comp
```

**Create The Shared Environment**

From the repository root:

```bash
conda env create -f environment.yml
conda activate tapemeasure-analyze
```

This creates one shared environment named `tapemeasure-analyze` containing everything needed for both `tapeseq` and `tapebulk`.

The environment file also installs this repository in editable mode with `pip install -e .`, so changes in the repository are reflected in the installed commands without needing to reinstall from scratch every time.

**Verify The Install**

```bash
tapeseq --help
tapebulk --help
```

If both commands print help text, the installation is working.

## Updating An Existing Local Install

If you already installed an older version of this repository, update the shared environment with:

```bash
conda env update -f environment.yml --prune
conda activate tapemeasure-analyze
pip install -e .
hash -r
```

Notes:

- `conda env update ... --prune` adds new dependencies and removes old ones not listed anymore
- `pip install -e .` refreshes the installed console entry points, including `tapebulk`
- `hash -r` clears the shell command cache in `bash` so the terminal notices any updated command paths

## How Output Is Organized

Both tools write results into an output directory you choose.

Each run writes shared metadata under:

- `OUTPUT_DIR/run_metadata/analysis_parameters.json`
- `OUTPUT_DIR/run_metadata/run.log`
- `OUTPUT_DIR/run_metadata/current_state.json`

Analysis-specific outputs are written under `per_accession/` or `global_outputs/`, depending on the tool and mode.

## `tapeseq`

`tapeseq` is the sequence-analysis side of the project.

### Core Idea

For each input sequence, `tapeseq` builds a self-similarity matrix using BLOSUM45 residue scores. It then smooths that matrix, converts it to a 1D signal by taking the row-wise maximum away from the diagonal, and standardizes that signal against shuffled versions of the same sequence. That standardized projection is the main signal reused by downstream sequence analyses.

Full-sequence artifacts are cached by default in `.tapeseq_cache/` so repeated runs do not need to rebuild the same projections.

### General `tapeseq` Usage

```bash
tapeseq MODE -i INPUT.fasta -o OUTPUT_DIR [options]
```

Main modes:

- `exploratory`
- `periodicity`
- `repeat-analysis`
- `map-domain`

You can inspect mode-specific arguments with:

```bash
tapeseq exploratory --help
tapeseq periodicity --help
tapeseq repeat-analysis --help
tapeseq map-domain --help
```

### `tapeseq exploratory`

Purpose:
Write diagnostic outputs showing how each sequence behaves under the self-similarity projection workflow.

What it does:

1. Loads sequences from FASTA.
2. Builds the full self-similarity matrix.
3. Smooths the matrix iteratively.
4. Computes the standardized projection from shuffled-sequence baselines.
5. Writes heatmaps and projection plots.

Typical command:

```bash
tapeseq exploratory -i INPUT.fasta -o OUTPUT_DIR
```

Optional trimmed analysis:

```bash
tapeseq exploratory -i INPUT.fasta -o OUTPUT_DIR --trimmings TRIMS.csv
```

If trimming rules are provided, the tool writes outputs for the trimmed windows instead of the full sequences.

Typical outputs:

- per-accession raw matrix heatmaps
- per-accession smoothed matrix heatmaps
- per-accession standardized projection plots
- `global_outputs/trimmed_sequences.fasta` when trims are used

### `tapeseq periodicity`

Purpose:
Estimate whether the standardized sequence projection contains a dominant repeating period.

What it does:

1. Reuses the standardized projection from the full sequence or trimmed region.
2. Runs an FFT on that 1D signal.
3. Searches a configurable frequency window for the strongest peak.
4. Converts that peak into an inferred period length.
5. Slices the analyzed sequence into putative repeat segments using that inferred period.

Typical command:

```bash
tapeseq periodicity -i INPUT.fasta -o OUTPUT_DIR --trimmings TRIMS.csv
```

Important interpretation note:
The repeat segments written by this mode are putative repeats derived from FFT-estimated period lengths. They are meant as hypothesis-generating outputs, not validated biological repeat annotations.

Typical outputs:

- per-accession projection plots
- per-accession FFT power plots
- per-accession FASTA files of putative segments
- `global_outputs/periodicity/results_summary.csv`
- `global_outputs/periodicity/putative_repeats_all.fasta`

### `tapeseq repeat-analysis`

Purpose:
Score externally supplied repeats against the standardized projection of the full sequence.

What it does:

1. Loads full sequences from FASTA.
2. Loads supplied repeat records from a repeat FASTA.
3. Locates each repeat in the full sequence by exact match or local alignment fallback.
4. Scores the located repeat window using the maximum value in the standardized projection over that interval.

Typical command:

```bash
tapeseq repeat-analysis -i INPUT.fasta -o OUTPUT_DIR --repeats-fasta REPEATS.fasta
```

Repeat FASTA headers should follow this pattern:

```text
ACCESSION_R1
ACCESSION_R2
ACCESSION_RF
```

Typical outputs:

- per-accession repeat window score tables
- `global_outputs/repeat_analysis/supplied_repeat_scores_summary.csv`

### `tapeseq map-domain`

Purpose:
Project a known domain from one aligned reference accession onto other aligned accessions and convert that mapping into conservative trimming rules.

What it does:

1. Loads an aligned FASTA.
2. Locates the reference domain inside the ungapped reference accession.
3. Converts that reference domain into alignment-column coordinates.
4. Projects those columns onto every accession in the alignment.
5. Emits accession-specific trimming rules only when the projected mapping is conservative enough.
6. Records QC for all accessions, including partial or weak mappings.

Typical command:

```bash
tapeseq map-domain \
  -a ALIGNMENT.fasta \
  -o OUTPUT_DIR \
  --reference-accession REF_ACC \
  --reference-domain DOMAIN.fasta
```

Typical outputs:

- `global_outputs/map_domain/mapped_domain_trimmings.csv`
- `global_outputs/map_domain/mapped_domain_qc.csv`

### Trimming CSV Format

For `exploratory` and `periodicity`, trimming rules must contain `start` and `end`, plus either `tmp_length` or `accession`.

Length-based example:

```csv
tmp_length,start,end
120,10,95
150,20,130
```

Accession-based example:

```csv
accession,start,end
ACC001,10,95
ACC002,20,130
```

If both are present, accession-specific trimming rules take precedence.

## `tapebulk`

`tapebulk` is the coordinate-analysis side of the project.

### Core Idea

Instead of looking at sequence self-similarity, `tapebulk` looks at where residues sit in 3D space. The current implementation expects a CSV containing one row per residue coordinate:

```text
Pos, Res, X, Y, Z
```

This can be used to ask whether some parts of a chain are locally crowded, locally expanded, or globally unusual in their pairwise distance distribution.

### General `tapebulk` Usage

```bash
tapebulk MODE -i POINTS.csv -o OUTPUT_DIR [options]
```

Main modes:

- `bulkiness`
- `pairwise-histogram`

Inspect arguments with:

```bash
tapebulk bulkiness --help
tapebulk pairwise-histogram --help
```

### `tapebulk bulkiness`

Purpose:
Measure local spatial crowding along the chain.

What it does:

1. Loads residue coordinates from CSV.
2. Uses a KD-tree for efficient neighborhood queries.
3. Counts neighbors within one or more radii for every residue.
4. Optionally subtracts the expected neighbor count for a perfectly straight chain using a configurable residue spacing.
5. Computes the mean distance to the first `k` nearest neighbors.
6. Converts that mean nearest-neighbor distance into a simple crowding score:

```text
bulkiness score = 1 / (mean k-nearest-neighbor distance)^3
```

Typical command:

```bash
tapebulk bulkiness \
  -i POINTS.csv \
  -o OUTPUT_DIR \
  --radii 6 10 15 \
  --k-neighbors 10 \
  --subtract-baseline \
  --make-plots
```

Typical outputs:

- `global_outputs/bulkiness/bulkiness_table.csv`
- `global_outputs/bulkiness/bulkiness_summary.json`
- optional PNG plots for neighbor counts, extra crowding, and kNN bulkiness

### `tapebulk pairwise-histogram`

Purpose:
Summarize the global distance distribution of the coordinate set.

What it does:

1. Loads residue coordinates from CSV.
2. Computes all pairwise Euclidean distances.
3. Writes the distance list to CSV.
4. Writes a histogram plot.
5. Writes summary statistics as JSON.

Typical command:

```bash
tapebulk pairwise-histogram -i POINTS.csv -o OUTPUT_DIR --histogram-bins 40
```

Typical outputs:

- `global_outputs/bulkiness/pairwise_distances.csv`
- `global_outputs/bulkiness/pairwise_distance_summary.json`
- `global_outputs/bulkiness/pairwise_distance_hist.png`

## Compatibility Wrappers

Older repository-local entrypoints still exist:

```bash
python3 unified_sequence_analysis.py --help
python3 bulkiness_v2.py --help
```

These wrappers forward into the packaged implementations and are mainly there to avoid breaking older local workflows.
