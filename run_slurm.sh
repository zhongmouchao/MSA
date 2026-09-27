#!/bin/bash
# =====================================================================
# SLURM submission script for the MSA + motif extraction + counting pipeline
#
# Usage:
#   cd to the pipeline folder (the one containing run_pipeline.py), then:
#       sbatch run_slurm.sh                    # default yearly counts
#       sbatch run_slurm.sh --monthly          # extra pipeline args pass through
#       sbatch run_slurm.sh --start-year 2000 --motif-min-called 1.0
#
# The pipeline processes every FASTA dataset in ./input (each with a
# matching .txt config, or the shared input/reference.txt) and writes
# results to ./output.
# =====================================================================
#SBATCH --job-name=msa_pipeline
#SBATCH --output=logs/msa_%j.out
#SBATCH --error=logs/msa_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4        # hmmalign is single-threaded; 4 is plenty for I/O overlap
#SBATCH --mem=16G                # ~550k sequences need ~6 GB; raise for larger inputs
#SBATCH --time=24:00:00
#SBATCH --partition=default_partition
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zc83@cornell.edu

set -euo pipefail

# --- activate the conda environment --------------------------------------
# Auto-detect the conda base from the `conda` command itself (works when the
# env lives in ~/.conda/envs rather than under the base install). Override
# with CONDA_BASE=/path at submission if conda is not on PATH on compute
# nodes, or replace these lines with:  module load anaconda && source activate msa
if [ -z "${CONDA_BASE:-}" ]; then
    if command -v conda >/dev/null 2>&1; then
        CONDA_BASE="$(conda info --base)"
    else
        CONDA_BASE="$HOME/miniconda3"
    fi
fi
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate msa

# --- run from the directory where sbatch was invoked ----------------------
cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

echo "Job $SLURM_JOB_ID on $SLURM_NODELIST started $(date)"
echo "Pipeline args: $*"

python3 run_pipeline.py "$@"

echo "Finished $(date)"
