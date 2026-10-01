# MSA + Motif Extraction + Counting Pipelines

This workspace hosts **two related workflows** that share the same core
engine (`run_pipeline.py`, untouched):

| | Workflow A: original batch pipeline | Workflow B: viroGym benchmark pipeline |
|---|---|---|
| Purpose | align a batch of FASTAs against ONE reference, extract motifs, count | produce DMS training data for the 79 viroGym benchmarks |
| Driver | `submit_list.sh` (+ `run_slurm.sh`) | `main_run.py` (chains 6 modules) |
| Input | flat `input/` folder | `RAW_FASTA_DATA/` (recursive) + `Virogym_Benchmark/` |
| Config | `input/reference.txt` or per-file `<stem>.txt` | auto-generated per group in `planned_configs/` |
| Output | `output/<chunk>/` | `aligned/<group>/<chunk>/` + `training_data/` |
| Job names | `msa_<chunk>` | `vg_<group>_<chunk>` |

The two systems never collide: different job-name prefixes and different
output roots, so both may run on the cluster at the same time.

---

## Workflow A — original batch pipeline

```
input/
  <chunk>.fasta            # surveillance sequences; header must contain a year
  reference.txt            # shared config (or <chunk>.txt per dataset)
output/
submit_list.sh             # bash submit_list.sh   (one SLURM job per FASTA)
run_slurm.sh               # per-job wrapper (4 CPU / 48G / 24h)
run_pipeline.py            # the engine; do not edit casually
```

Config file format (`reference.txt`):

```
>protein_reference
<whole-protein amino-acid sequence>          # required
>motif_reference_<tag>
<motif amino-acid sequence>                  # optional, repeatable
```

- protein reference only  -> align + count, no motif extraction
- one or more motifs      -> align once, extract + count each motif separately
- motif identical to the protein reference ("whole") -> count only sequences
  with <=10% gap/X across the whole protein (`--motif-min-called 0.9` default)

Manual single job:  `sbatch run_slurm.sh --dataset <chunk.fasta> --monthly`

Useful options: `--monthly` (YYYY-MM columns), `--start-year N`,
`--motif-min-called 0..1`, `--min-orf-frac`, `--min-orf-aa`.

Resume: `submit_list.sh` skips chunks with `output/<chunk>/.msa_complete`
and chunks already in the queue. Splitting oversized inputs first:
see `split_fasta.py` / `submit_split.sh` (if present).

---

## Workflow B — viroGym benchmark pipeline

```
Virogym_Benchmark/
  virogym_79_benchmarks_combined_with_target_type.csv   # the registry (inside!)
  <79 raw DMS files>.csv
RAW_FASTA_DATA/
  SCV2/                       # flat: genome + spike chunk fastas
  IAV/human/HA/  IAV/avian/HA/  IAV/avian/PB2/  IAV/all/...
  <VIRUS_CODE>/               # one folder per other virus, when you get data
references/                   # only for full_protein_plus_motif benchmarks
aligned/                      # created automatically; one folder per group
planned_configs/              # auto-generated run_pipeline configs
training_data/                # final per-benchmark training tables
check_registry.py scan_inputs.py scan_aligned.py plan_tasks.py
run_tasks.py aggregate_counts.py main_run.py
vg_run_slurm.sh               # per-job wrapper (4 CPU / 48G / 24h)
run_pipeline.py               # shared engine (same file as Workflow A)
```

### Daily use

```bash
python3 main_run.py --dry-run     # validate + plan + show wave plan (nothing submits)
python3 main_run.py               # full run: submit SLURM waves, wait, aggregate
python3 main_run.py --no-submit   # skip submission; re-aggregate finished chunks
```

Or step by step:

```bash
python3 check_registry.py     # step 0: folder <-> registry agreement
python3 scan_inputs.py        # step 1: -> input_inventory.json
python3 scan_aligned.py       # step 1.5: -> aligned_inventory.json
python3 plan_tasks.py         # step 2: -> task_plan.json + planned_configs/
python3 run_tasks.py          # step 3: submit + babysit SLURM waves
python3 aggregate_counts.py   # step 4: -> training_data/
```

### How it decides what to do

1. **Registry check** — every DMS file needs a registry row (with
   `target_seq`, `host`, `virus_code`, `recommended_reference`); new files
   are reported as `[NEW]` and skipped until you add their row.
2. **Grouping** — benchmarks sharing an identical de-gapped target sequence
   are aligned ONCE; different motifs on the same reference share the
   alignment and are counted separately.
3. **File suitability** — a FASTA serves a reference of length L aa if it is
   gene-length-matched (0.6–1.5 x L), parent-protein-matched (for region
   targets like RBD inside spike), or genome-length (>=8000 nt).
   IAV routing: most-specific host folder wins (`human/`,`avian/` before
   `all/`, never mixed); filename segment hints (H1/H3/PB2/...) are enforced.
4. **Wave scheduler (cluster politeness)** — one SLURM job per
   (group x FASTA). Multi-file groups run EXCLUSIVELY (all their jobs in
   parallel, nothing else); single-file groups run up to
   `--parallel-small` (default 10) at a time.
5. **Aggregation** — per-benchmark `training_data/<benchmark>_training_counts.csv`.
   Dates are strictly validated (SCV2 >= 2019, HIV >= 1980, flu >= 1900,
   others >= 1950, <= next year); anything else folds into ONE `unknown`
   column. Variants with >10% `-`/`X` in the motif peptide are dropped.
   Audit in `training_data/_aggregation_summary.csv`.

### Resume / incremental behavior

- Everything is state on disk: `aligned/<group>/<chunk>/.msa_complete`.
- Interrupted run? Re-run the same command — finished chunks are skipped.
- Added a DMS file? Drop it in `Virogym_Benchmark/`, add its registry row,
  re-run `main_run.py` — only the new work executes.
- Added FASTAs? Same — only the new chunks align.
- Changed your mind about a group? Delete its `aligned/<group>/` folder and
  re-run; only that group redoes.

### Pre-existing alignments (optional)

```bash
python3 scan_aligned.py --register aligned/<old_folder> --config <old reference.txt>
```

Registers an old alignment for reuse ONLY if its reference is identical
(hash-matched) to a group reference; also reports motif-level matches whose
counts can be reused directly. Reuse is exact-match only — never coordinate
mapping.

### Special cases

- `PESV_stability` is the only benchmark (1/79) needing
  `full_protein_plus_motif`: supply the full-length Q9QEJ5 polyprotein as
  `references/PESV_stability_fullprotein.fasta`, else it stays blocked.
- SCV2 benchmarks bin by month; all other viruses by year.

### Cluster environment

```bash
conda env create -f environment.yml   # once: creates env 'msa'
conda activate msa                    # python3, pandas, biopython, hmmer
```

---

## Troubleshooting

| Symptom | Action |
|---|---|
| `sbatch: ... Socket timed out` | scheduler hiccup; re-run the same command later — submitted jobs are detected by name and not duplicated |
| chunk job finished but no `.msa_complete` | read `logs/vg_<jobid>.out/.err` and the chunk `.log`; fix input, re-run `run_tasks.py` (only missing chunks resubmit) |
| benchmark skipped "no FASTA folder matched" | add `RAW_FASTA_DATA/<virus_code>/...` data, re-run |
| benchmark blocked (PESV) | provide `references/PESV_stability_fullprotein.fasta` |
| memory failure on huge genome jobs | split the FASTA (Workflow A split scripts), or raise `--mem` in `vg_run_slurm.sh` |
