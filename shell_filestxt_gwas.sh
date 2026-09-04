#!/bin/bash
#$ -cwd
#$ -N GWAS_filestxt
#$ -l h_rt=24:00:00
#$ -l h_data=12G
#$ -pe shared 8
#$ -o /u/scratch/v/vturnbil/GWAS_preprocessed/logs/joblog.$JOB_ID
#$ -e /u/scratch/v/vturnbil/GWAS_preprocessed/logs/joblog.$JOB_ID
#$ -M vturnbil@g.ucla.edu
#$ -m bea

echo "Job $JOB_ID started on $(hostname) at $(date)"

# -------------------------
# Load conda properly
# -------------------------
source /u/project/xyang123/vturnbil/packages/miniconda3/etc/profile.d/conda.sh

# Activate scing environment
conda activate scing

# Safety check (optional but recommended)
echo "Using python:"
which python
python -c "import sys; print(sys.executable)"

# -------------------------
# Run job
# -------------------------
cd /u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MDF || exit 1

python filestxt_gwas.py

echo "Job $JOB_ID finished at $(date)"
