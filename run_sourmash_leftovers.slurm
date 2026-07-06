#!/bin/bash
#SBATCH -J sourmash_leftovers_MW153array
#SBATCH --account=gts-mwhiteley3-biocluster
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o sourmash_leftovers_MW153array-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hwilms3@gatech.edu

set -euo pipefail

# Keep any numerical libraries from trying to use extra cores per process.
# The Python script itself controls sample-level parallelism.
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
echo "Sourmash: $(which sourmash)"
echo "SLURM CPUs: ${SLURM_CPUS_PER_TASK}"
echo "========================================"

python run_sourmash_on_non_pao1_reads.py \
  --config sourmash_config.yaml

echo "========================================"
echo "Job finished: $(date)"
echo "========================================"
