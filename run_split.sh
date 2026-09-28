#!/bin/bash
# =====================================================================
# SLURM job: split ONE large FASTA into N parts (phase 1).
# Called by submit_list.sh as:  sbatch run_split.sh "<fasta file>" [parts]
# =====================================================================
#SBATCH --job-name=msa_split
#SBATCH --output=logs/split_%j.out
#SBATCH --error=logs/split_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --partition=default_partition
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=zc83@cornell.edu

set -euo pipefail

if [ -z "${CONDA_BASE:-}" ]; then
    if command -v conda >/dev/null 2>&1; then
        CONDA_BASE="$(conda info --base)"
    else
        CONDA_BASE="$HOME/miniconda3"
    fi
fi
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate msa

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

FASTA="$1"
PARTS="${2:-4}"

echo "Job $SLURM_JOB_ID on $SLURM_NODELIST started $(date)"
echo "Splitting input/$FASTA into $PARTS parts"

python3 split_fasta.py "input/$FASTA" --parts "$PARTS" --out-dir input

echo "Finished $(date)"
