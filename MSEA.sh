#!/bin/bash
#$ -cwd
#$ -N MSEA
#$ -l highp
#$ -l h_rt=72:00:00
#$ -l h_data=32G
#$ -o /u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MSEA/logs/joblog.$JOB_ID.$TASK_ID.out
#$ -e /u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MSEA/logs/joblog.$JOB_ID.$TASK_ID.err
#$ -M vturnbil@g.ucla.edu
#$ -m bea
#$ -j y
#$ -t 1-202

source /u/local/Modules/default/init/modules.sh
module load apptainer

# ---- get list of GWAS folders automatically ----
TRAITS=(/u/scratch/v/vturnbil/GWAS_post_mdf/*)

TRAIT=$(basename "${TRAITS[$SGE_TASK_ID-1]}")

echo "Running trait $TRAIT"
date
hostname

apptainer exec \
/u/local/apps/apptainer/H2_containers/h2-rstudio_4.5.0.sif \
Rscript \
/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MSEA/run_msea.R \
"$TRAIT"

echo "Finished $TRAIT"
date
