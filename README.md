# Multiple Sequence Alignment Pipeline

Streamlined, cluster-ready version of the original notebook
(`26zc0101_Sequence Alignment_with known nucleotides.ipynb`).
One runner loops over every dataset in `input/`, auto-detects what kind of
data it is, runs the right stages, and writes results to `output/`.

## Requirements

- Python 3.9+ with `biopython` and `pandas`
- HMMER (`hmmbuild`, `hmmalign`) on `PATH` — e.g. `brew install hmmer` or
  `module load hmmer` on a cluster

## How to run

```bash
cd "Multiple Sequence Alignment"
python3 run_pipeline.py                  # process every dataset in ./input
python3 run_pipeline.py --monthly        # split counts by year AND month
python3 run_pipeline.py --input-dir mydata --output-dir myresults
```

Options: `--monthly`, `--start-year N`, `--min-orf-frac F` (minimum ORF/CDS
length as a fraction of the reference length, default 0.5) or `--min-orf-aa N`
(absolute override), `--input-dir`, `--output-dir`.

On a cluster, submit it as a normal job, e.g.:

```bash
#!/bin/bash
#SBATCH --mem=32G --time=12:00:00
cd /path/to/"Multiple Sequence Alignment"
python3 run_pipeline.py
```

## Preparing a dataset

Put **two files with the same base name** into `input/`:

```
input/
    my_virus.fasta     # your sequences
    my_virus.txt       # config: reference sequence(s), FASTA-style
```

If many datasets share the same references, you can instead drop a single
`reference.txt` (or `config.txt`) into `input/` — any dataset without its own
`<name>.txt` will use the shared one.

Config file (`my_virus.txt`):

```
>protein_reference
MFVFLVLLPLVSSQCVNLTTRTQ...
>motif_reference
DPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKF
```

What you put in the config tells the pipeline what you want:

| Config contains       | What runs                                            |
|-----------------------|------------------------------------------------------|
| protein + motif       | ORF (if genomes) → MSA → motif extraction → counting |
| protein only          | ORF (if genomes) → MSA → counting (full sequences)   |
| motif only            | motif extraction → counting (input must be aligned)  |

Sequence headers must contain a collection date somewhere:
`YYYY`, `YYYY-MM`, or `YYYY-MM-DD` (e.g. GISAID `>hCoV-19/...|2023-05-12`).
Records without a date land in the `unknown` column.

## What the pipeline decides automatically

For each dataset it samples the first 1,000 sequences and detects:

- **nucleotide vs protein** — from the alphabet
- **raw vs already aligned** — gaps present + uniform length = aligned
  (aligned input skips the MSA entirely)
- **whole genome vs gene/CDS** — median length vs. the protein reference
  (genomes get ORF extraction first: the ORF with the highest k-mer coverage
  of the reference is picked, so the reference may be a full-length protein
  or just a fragment of it; the minimum ORF length defaults to 50% of the
  reference length — `--min-orf-frac` / `--min-orf-aa`)

It also cross-checks your config against the data and stops with a clear
message on impossible combinations (e.g. motif-only config + raw
unaligned input, which cannot be aligned without the protein reference).

## Outputs (per dataset, in `output/<name>/`)

| File | Content |
|---|---|
| `<name>_orf_nucleotides.fasta` / `_orf_proteins.fasta` | ORF stage (genome input only) |
| `<name>_aligned_proteins.fasta` | ref-only HMM protein alignment (incl. reference row) |
| `<name>_aligned_nucleotides.fasta` | codon-aware back-translated NT alignment |
| `<name>_motif_proteins.fasta` / `_motif_nucleotides.fasta` | paired motif fragments (if motif config given) |
| `<name>_counts_by_year_wide.csv` | always produced; `--monthly` gives `_by_month_wide.csv` |
| `<name>.log` | full log of that dataset's run |

Counts CSV format: `Motif Nucleotide, Motif Peptide, before YYYY, YYYY, …,
unknown, Count` — one row per unique variant. In `--monthly` mode the date
columns are `YYYY-MM` plus a `YYYY-unknown` bucket per year for records that
only carry a year.

## Files

- `run_pipeline.py` — the pipeline (this is all you need on the cluster)
- `input/` — put your FASTA + config pairs here (`example_demo.*` is a tiny
  working example you can delete)
- `output/` — results
- `make_test_data.py` — regenerates the synthetic example data
- the original notebook and an example counts CSV are kept here for reference
