#!/bin/bash
#SBATCH -J extract_non_decoy_non_pao1_reads_MW153array
#SBATCH --account=gts-mwhiteley3-biocluster
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -t 12:00:00
#SBATCH -o extract_non_decoy_non_pao1_reads_MW153array-%j.out
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

python extract_non_decoy_non_pao1_reads.py \
    --config extract_non_pao1_config.yaml

echo "========================================"
echo "Job finished: $(date)"
echo "========================================"