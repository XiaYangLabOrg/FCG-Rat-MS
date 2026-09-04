#!/bin/bash
#$ -cwd
#$ -j y
#$ -V

#$ -N SATURN
#$ -o /u/project/xyang123/xyang123-NOBACKUP/turnbill/Tools/SATURN/logs/SATURN.$JOB_ID

#$ -l L40S=1
#$ -l gpu=1
#$ -l h_data=24G
#$ -l h_rt=23:59:00

#$ -pe shared 1

#$ -M $USER@mail
#$ -m bea

source /u/local/apps/mambaforge/23.11.0/etc/profile.d/conda.sh
conda activate saturn_saturn

export NUMBA_CACHE_DIR=/tmp/$USER/numba_cache
mkdir -p $NUMBA_CACHE_DIR

echo "Host: $(hostname)"
echo "Date: $(date)"
echo "Python:"
which python
python --version

nvidia-smi

cd /u/project/xyang123/xyang123-NOBACKUP/turnbill/Tools/SATURN

python train-saturn.py \
    --in_data in_data.csv \
    --work_dir SATURN_output \
    --device_num 0 \
    --ref_label_col celltype
