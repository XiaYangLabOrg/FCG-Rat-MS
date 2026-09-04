#!/bin/bash
#$ -cwd
#$ -o /u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/logs_wKDA/joblog.$JOB_ID.$SGE_TASK_ID
#$ -j y
#$ -l h_rt=24:00:00,h_data=24G
#$ -pe shared 8
#$ -M $USER@mail
#$ -m bea
#$ -N wKDA
#$ -t 1-8

echo "======================================="
echo "Job ID: $JOB_ID"
echo "Task ID: $SGE_TASK_ID"
echo "Host: $(hostname)"
date
echo "======================================="

# Initialize module system
source /u/local/Modules/default/init/modules.sh

# Clean environment
module purge
module load apptainer

# ---------------------------------------
# LIST OF NETWORK DIRECTORIES
# ---------------------------------------

DATASETS=(
/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/MSMC_Astrocytes
/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/MSMC_Cholinergic_Neurons
/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/MSMC_GABAergic_Neurons_1
/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/MSMC_GABAergic_Neurons_2
/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/MSMC_Glutamatergic_Neurons
/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/MSMC_Immature_Oligodendrocytes
/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/MSMC_Mature_Oligodendrocytes
/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output/MSMC_OPCs
)

# Shared DEG file
DEG_FILE="/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Limma/deg_all2.RDS"

# ---------------------------------------
# SELECT DATASET FOR THIS TASK
# ---------------------------------------

BASE_DIR=${DATASETS[$((SGE_TASK_ID-1))]}

if [ -z "$BASE_DIR" ]; then
    echo "ERROR: No dataset defined for task $SGE_TASK_ID"
    exit 1
fi

echo "Running dataset:"
echo "$BASE_DIR"
echo "DEG file:"
echo "$DEG_FILE"
echo "---------------------------------------"

# ---------------------------------------
# RUN R SCRIPT INSIDE CONTAINER
# ---------------------------------------

apptainer exec \
/u/local/apps/apptainer/H2_containers/h2-rstudio_4.5.0.sif \
Rscript \
/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/R_Files/wKDA_after_SCING.R \
"$BASE_DIR" "$DEG_FILE"

echo "Finished task $SGE_TASK_ID"
date
echo "======================================="
