#!/bin/bash
#$ -cwd
#$ -j y
#$ -V
#$ -o /u/project/xyang123/xyang123-NOBACKUP/turnbill/Tools/SATURN/logs/mouse_esm.$JOB_ID
#$ -N mouse_esm
#$ -M $USER@mail
#$ -m bea

#$ -l h_data=24G
#$ -pe shared 4
#$ -l h_rt=23:59:00

source /u/local/apps/mambaforge/23.11.0/etc/profile.d/conda.sh
conda activate saturn_saturn

cd /u/project/xyang123/xyang123-NOBACKUP/turnbill/Tools/SATURN/esm

python ./scripts/extract.py esm1b_t33_650M_UR50S \
/u/project/xyang123/xyang123-NOBACKUP/turnbill/Tools/SATURN/Proteins/Mus_musculus.GRCm39.pep.all_clean.fa \
/u/scratch/v/vturnbil/SATURN_esm/mouse \
--include mean --truncate
