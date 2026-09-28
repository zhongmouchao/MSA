#!/bin/bash
# Run one FASTA through the existing MSA + motif extraction + counting pipeline.
# Normally submitted by: bash submit_list.sh
# Individual job (create logs/ first):
#   sbatch run_slurm.sh --dataset <chunk.fasta> --monthly
#SBATCH --job-name=msa_pipeline
#SBATCH --output=logs/msa_%j.out
#SBATCH --error=logs/msa_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --partition=default_partition
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zc83@cornell.edu

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

if [ "$#" -lt 2 ] || [ "$1" != "--dataset" ]; then
    echo "Usage: sbatch run_slurm.sh --dataset <chunk.fasta> [pipeline options]" >&2
    exit 1
fi
FASTA="$2"
shift 2
if [[ "$FASTA" == */* ]] || [ ! -s "input/$FASTA" ]; then
    echo "Expected the filename of a nonempty FASTA in input/: $FASTA" >&2
    exit 1
fi

# Keep job-specific input/output selection fixed; other pipeline options pass through.
for arg in "$@"; do
    case "$arg" in
        --dataset|--dataset=*|--input-dir|--input-dir=*|--output-dir|--output-dir=*)
            echo "Do not override --dataset, --input-dir, or --output-dir here." >&2
            exit 1 ;;
    esac
done

BASE="${FASTA%.*}"
JOB_OUTPUT="$SLURM_SUBMIT_DIR/output/$BASE"
mkdir -p "$JOB_OUTPUT"

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

echo "Job $SLURM_JOB_ID on $SLURM_NODELIST started $(date)"
echo "Dataset: $FASTA"
echo "Job output folder: $JOB_OUTPUT"

# The Python pipeline adds a dataset subfolder under --output-dir.
# A separate output root per chunk prevents concurrent jobs from overwriting
# run_summary.csv and integrated count tables. These tables cover this chunk only.
rm -f "$JOB_OUTPUT/.msa_complete"
python3 run_pipeline.py \
    --dataset "$FASTA" \
    --input-dir "$SLURM_SUBMIT_DIR/input" \
    --output-dir "$JOB_OUTPUT" \
    "$@"

touch "$JOB_OUTPUT/.msa_complete"
echo "Finished $(date)"
