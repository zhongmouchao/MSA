#!/bin/bash
# =====================================================================
# vg_run_slurm.sh -- SLURM wrapper for ONE (reference group x FASTA) job
# of the viroGym surveillance pipeline.
#
# Submitted by run_tasks.py; normally not run by hand.
# Manual use (from the pipeline root, after mkdir -p logs):
#   sbatch vg_run_slurm.sh --input-dir jobs/<g>/<chunk>/input \
#       --output-dir aligned/<g>/<chunk> --dataset <file.fasta> --monthly
# =====================================================================
#SBATCH --job-name=vg_pipeline
#SBATCH --output=logs/vg_%j.out
#SBATCH --error=logs/vg_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --partition=default_partition
## SBATCH --mail-type=ALL
## SBATCH --mail-user=zc83@cornell.edu

set -euo pipefail

INPUT_DIR=""
OUTPUT_DIR=""
DATASET=""
EXTRA=()
while [ "$#" -gt 0 ]; do
    case "$1" in
        --input-dir)  INPUT_DIR="$2";  shift 2 ;;
        --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
        --dataset)    DATASET="$2";    shift 2 ;;
        *)            EXTRA+=("$1");   shift ;;
    esac
done
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$DATASET" ]; then
    echo "Usage: vg_run_slurm.sh --input-dir D --output-dir D --dataset F [pipeline opts]" >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"
mkdir -p logs "$OUTPUT_DIR"

# Prefer cluster scratch; /tmp is a fallback and may be RAM-backed.
export TMPDIR="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
mkdir -p "$TMPDIR"

if [ -z "${CONDA_BASE:-}" ]; then
    if command -v conda >/dev/null 2>&1; then
        CONDA_BASE="$(conda info --base)"
    else
        CONDA_BASE="$HOME/miniconda3"
    fi
fi
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate msa

echo "Job ${SLURM_JOB_ID:-manual} on ${SLURM_NODELIST:-localhost} started $(date)"
echo "Dataset: $DATASET"
echo "Input dir: $INPUT_DIR"
echo "Output dir: $OUTPUT_DIR"
echo "Pipeline args: ${EXTRA[*]:-}"

# The completion marker is removed first: it only exists if the WHOLE
# pipeline (align + extract + count) finished for this chunk.
rm -f "$OUTPUT_DIR/.msa_complete"
python3 "$SCRIPT_DIR/run_pipeline.py" \
    --dataset "$DATASET" \
    --input-dir "$INPUT_DIR" \
    --output-dir "$OUTPUT_DIR" \
    ${EXTRA[@]+"${EXTRA[@]}"}

touch "$OUTPUT_DIR/.msa_complete"
echo "Finished $(date)"
