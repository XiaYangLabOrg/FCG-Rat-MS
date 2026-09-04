#!/bin/bash
#$ -cwd
#$ -o /u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/logs/joblog.$JOB_ID.$SGE_TASK_ID
#$ -j y
#$ -V
#$ -l h_rt=47:00:00,h_data=24G,highp
#$ -pe shared 8
#$ -M $USER@mail
#$ -m bea
#$ -N SCING
#$ -t 1-24

LOGDIR=/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/logs
WORKDIR=/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/adata_by_celltype
FILELIST=/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/filelist.txt

mkdir -p "$LOGDIR"
cd "$WORKDIR" || exit 1

echo "Job $JOB_ID Task $SGE_TASK_ID started on $(hostname) at $(date)"

# -------------------------
# Conda
# -------------------------
source /u/project/xyang123/vturnbil/packages/miniconda3/etc/profile.d/conda.sh
conda activate scing

# -------------------------
# Get input file for this task
# -------------------------
FILE=$(sed -n "${SGE_TASK_ID}p" "$FILELIST")

if [ -z "$FILE" ]; then
    echo "ERROR: No file found for task $SGE_TASK_ID in $FILELIST"
    exit 1
fi

# Optional: if filenames are relative, convert to full path
# FILE="$WORKDIR/$FILE"

if [ ! -f "$FILE" ]; then
    echo "ERROR: File does not exist: $FILE"
    exit 1
fi

echo "Processing file: $FILE"

# -------------------------
# Run SCING
# -------------------------
python "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/run_scing_MSMC.py" "$FILE"
STATUS=$?

echo "Job $JOB_ID Task $SGE_TASK_ID finished at $(date) with status $STATUS"
exit $STATUS
