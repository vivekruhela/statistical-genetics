#!/usr/bin/env bash
#SBATCH --job-name=statgen-build
#SBATCH --cpus-per-task=4
#SBATCH --mem=14G
#SBATCH --time=12:00:00
#SBATCH --output=statgen-build-%j.log
# The full data build (and, with RUN_TESTS=1, the checks) as a Slurm job on a compute node.
# Submit from the repository root:
#   sbatch simulate/run_all.sbatch data                    # the whole build
#   FROM=05 sbatch simulate/run_all.sbatch data            # resume at step 5, keep steps 0-4
#   FROM=05 RUN_TESTS=1 sbatch simulate/run_all.sbatch data
# Follow it with: tail -f statgen-build-<jobid>.log
# Edit the #SBATCH lines to suit your cluster (add --partition=... if it needs one).
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate statgen-r
echo "running on $(hostname), $(date)"
THREADS=${SLURM_CPUS_PER_TASK:-2} bash simulate/run_all.sh "${1:-data}"
if [[ ${RUN_TESTS:-0} == 1 ]]; then bash tests/run_tests.sh "${1:-data}"; fi
