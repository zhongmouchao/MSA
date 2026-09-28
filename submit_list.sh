#!/bin/bash
# =====================================================================
# Submit one SLURM job per FASTA dataset, running ONE AT A TIME.
#
# 1. Edit the DATASETS list below: one FASTA file name per line
#    (files must be in ./input, together with reference.txt or a
#    per-dataset <name>.txt config).
# 2. From the pipeline folder:   bash submit_list.sh
#
# Each item becomes:  sbatch run_slurm.sh --dataset "<file>"
# --dependency=singleton (same job name) makes the cluster run them
# sequentially. After every job finishes, the pipeline rebuilds the
# integrated count tables in ./output from everything completed so far.
# =====================================================================
set -euo pipefail
cd "$(dirname "$0")"

DATASETS=(
    "spikenuc0719_part_1.fasta"
    "spikenuc0719_part_2.fasta"
    "spikenuc0719_part_3.fasta"
    "spikenuc0719_part_4.fasta"
    "spikenuc0719_part_5.fasta"
    "spikenuc0719_part_6.fasta"
    "spikenuc0719_part_7.fasta"
    "spikenuc0719_part_8.fasta"
    "spikenuc0719_part_9.fasta"
    "spikenuc0719_part_10.fasta"
    "spikenuc0719_part_11.fasta"
    "spikenuc0719_part_12.fasta"
    "spikenuc0719_part_13.fasta"
    "spikenuc0719_part_14.fasta"
    "spikenuc0719_part_15.fasta"
    "spikenuc0719_part_16.fasta"
    "spikenuc0719_part_17.fasta"
    "spikenuc0719_part_18.fasta"

)

EXTRA_ARGS=(--monthly)          # e.g. EXTRA_ARGS=(--monthly --start-year 2000)

for f in "${DATASETS[@]}"; do
    if [ ! -f "input/$f" ]; then
        echo "SKIP: input/$f not found" >&2
        continue
    fi
    echo "Submitting: $f"
    sbatch run_slurm.sh --dataset "$f" "${EXTRA_ARGS[@]}"
done

echo "All ${#DATASETS[@]} dataset(s) submitted (singleton queue: one at a time)."
