#!/bin/bash
# Submit one independent MSA/motif/counting job per .fasta file in input/.
# Put this script and run_slurm.sh alongside run_pipeline.py in ~/MSA.
# Run: bash submit_list.sh
# No split jobs are submitted. Monthly counts are enabled below.
set -euo pipefail
cd "$(dirname "$0")"

if [ "$#" -ne 0 ]; then
    echo "Usage: bash submit_list.sh (no split/run argument needed)" >&2
    exit 1
fi

EXTRA_ARGS=(--monthly)
mkdir -p logs

# Query once. If Slurm cannot answer, do not mistake that for an empty queue.
if ! queued_jobs=$(squeue -u "$USER" -h -o '%j'); then
    echo "Could not check queued jobs. No jobs submitted; try again when Slurm responds." >&2
    exit 1
fi

found=0
submitted=0
skipped=0
failed=0

for chunk in input/*.fasta; do
    [ -f "$chunk" ] || continue
    found=$((found + 1))
    filename="${chunk##*/}"
    base="${filename%.fasta}"
    job_name="msa_$base"

    if [ ! -s "$chunk" ]; then
        echo "SKIP: $filename is empty."
        skipped=$((skipped + 1))
        continue
    fi

    if [ -f "output/$base/.msa_complete" ]; then
        echo "SKIP: $filename completed successfully on an earlier run."
        skipped=$((skipped + 1))
        continue
    fi

    if printf '%s\n' "$queued_jobs" | grep -Fxq "$job_name"; then
        echo "SKIP: $filename already has a queued/running job."
        skipped=$((skipped + 1))
        continue
    fi

    echo "Submitting MSA job: $filename"
    if sbatch --job-name="$job_name" run_slurm.sh --dataset "$filename" "${EXTRA_ARGS[@]}"; then
        submitted=$((submitted + 1))
    else
        echo "WARNING: Submission not confirmed for $filename; check Slurm before retrying." >&2
        failed=$((failed + 1))
    fi
    sleep 1
done

if [ "$found" -eq 0 ]; then
    echo "No FASTA files found: input/*.fasta" >&2
    exit 1
fi

echo "Found: $found; submissions confirmed: $submitted; skipped: $skipped; unconfirmed: $failed."
echo "Check jobs with: squeue -u \"\$USER\""
[ "$failed" -eq 0 ]
