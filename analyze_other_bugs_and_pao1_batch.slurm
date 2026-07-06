#!/bin/bash
#SBATCH -J PAO1_analysis
#SBATCH --account=gts-mwhiteley3-biocluster
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o PAO1_analysis-%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=hwilms3@gatech.edu

source /storage/home/hcoda1/2/hwilms3/miniforge3/etc/profile.d/conda.sh
conda activate bio

cd /storage/home/hcoda1/2/hwilms3/r-mwhiteley3-0/whiteley_PA_pipeline_post_analysis

echo "========================================"
echo "Job started: $(date)"
echo "Running on: $(hostname)"
echo "Python: $(which python)"
echo "========================================"

python analyze_other_bugs_and_pao1_per_sample.py \
    --config config.yaml

echo "========================================"
echo "Job finished: $(date)"
echo "========================================"