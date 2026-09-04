#!/bin/bash
#$ -l h_rt=24:00:00
#$ -l h_data=10G
#$ -o /u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MDF/logs/
#$ -j y
#$ -M $USER@mail
#$ -m bea

INPUTDIR="/u/scratch/v/vturnbil/GWAS_preprocessed"

for LOCFILE in "$INPUTDIR"/*.txt
do
    i=$(basename "$LOCFILE" .txt)
    echo "Processing $i"

    GENFILE="/u/project/xyang123/rainyliu/temp/gwas_mdf/hg38.dist50.txt"
    LNKFILE="/u/project/xyang123/shared/datasets/LD_Mergeomics_Ready/LD50.CEU.txt"

    OUTPATH="/u/scratch/v/vturnbil/GWAS_post_mdf/${i}"
    mkdir -p "$OUTPATH"

    NTOP=0.5
    echo -e "MARKER\tVALUE" > /tmp/header.txt
    sort -r -g -k 2 "$LOCFILE" > /tmp/sorted.txt
    NMARKER=$(wc -l < /tmp/sorted.txt)
    NMAX=$(echo "($NTOP*$NMARKER)/1" | bc)
    head -n "$NMAX" /tmp/sorted.txt > /tmp/top.txt
    cat /tmp/header.txt /tmp/top.txt > /tmp/subset.txt

    echo "Start MDF"

    nice /u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MDF/mdprune \
        /tmp/subset.txt "$GENFILE" "$LNKFILE" "$OUTPATH"

done
