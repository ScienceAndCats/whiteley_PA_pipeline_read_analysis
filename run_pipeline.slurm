#!/bin/bash
#SBATCH -J PA_read_filter_pipeline
#SBATCH --account=gts-mwhiteley3-biocluster
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o PA_read_filter_pipeline-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hwilms3@gatech.edu

set -euo pipefail

# Keep numerical libraries from using extra cores per process. The Python
# pipeline controls sample-level parallelism with settings.max_parallel_samples.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

source /storage/home/hcoda1/2/hwilms3/miniforge3/etc/profile.d/conda.sh
conda activate bio

cd /storage/home/hcoda1/2/hwilms3/r-mwhiteley3-0/whiteley_PA_pipeline_post_analysis

echo "========================================"
echo "Job started: $(date)"
echo "Running on: $(hostname)"
echo "Python: $(which python)"
echo "Sourmash: $(command -v sourmash || true)"
echo "SLURM CPUs: ${SLURM_CPUS_PER_TASK}"
echo "========================================"

python pa_read_filter_pipeline.py --config pipeline_config.yaml

echo "========================================"
echo "Job finished: $(date)"
echo "========================================"
