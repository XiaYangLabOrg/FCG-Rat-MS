# %%
import csv
import re
import os
import math
import sys
from operator import itemgetter

# %%
def read_gwas(gwas_dataset, filename, p_index, chr_index, bp_index):
    gwas_dic = {}

    # ---- delimiter-safe GWAS parsing ----
    with open(gwas_dataset) as gwas_file:
        for line in gwas_file:
            try:
                row = line.strip().split()

                if (row[p_index] != "NA"):
                    if float(row[p_index]) > 0 and float(row[p_index]) < 1:
                        gwas_dic['chr' + row[chr_index] + ':' + row[bp_index]] = \
                            -math.log10(float(row[p_index]))
            except:
                continue

    print(
        'Read in GWAS complete. A total of '
        + str(len(gwas_dic))
        + ' SNPs.'
    )

    snplist = os.listdir(
        '/u/project/xyang123/shared/datasets/mergeomics_gwas_update_2024/chr.snp.hg19/'
    )

    with open('/u/scratch/v/vturnbil/GWAS_preprocessed/' + filename, 'w') as output:
        writer = csv.writer(output, delimiter='\t')
        writer.writerow(['MARKER', 'VALUE'])

        for item in snplist:
            print(item)
            with open(
                '/u/project/xyang123/shared/datasets/mergeomics_gwas_update_2024/chr.snp.hg19/' + item
            ) as snp_file:
                snp = csv.reader(snp_file, delimiter='\t')
                for row in snp:
                    try:
                        if row[0] + ':' + row[1] in gwas_dic:
                            writer.writerow(
                                [row[2], gwas_dic[row[0] + ':' + row[1]]]
                            )
                    except:
                        continue


def parse_file_line(line):
    """
    Parse a line like:
    GWAS_NAME(1,2,5)
    """
    base, indices = line.split('(')
    indices = indices.rstrip(')')
    chr_index, bp_index, pvalue_index = map(int, indices.split(','))
    return base, chr_index, bp_index, pvalue_index


# %%
# ---- USER SETTINGS ----
gwas_dir = '/u/scratch/v/vturnbil/GWAS'
output_dir = '/u/scratch/v/vturnbil/GWAS_preprocessed'
file_list = '/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MDF/files.txt'

os.makedirs(output_dir, exist_ok=True)

# ---- LOOP OVER GWAS FILES ----
with open(file_list) as f:
    gwas_entries = [line.strip() for line in f if line.strip()]

for entry in gwas_entries:
    base, chr_index, bp_index, pvalue_index = parse_file_line(entry)

    dataset_input = os.path.join(gwas_dir, base + '.tsv')
    name = base + '.txt'

    if not os.path.exists(dataset_input):
        print('SKIPPING (file not found):', base)
        continue

    print(
        f'\nProcessing GWAS: {base} '
        f'(chr={chr_index}, bp={bp_index}, p={pvalue_index})'
    )

    read_gwas(
        dataset_input,
        name,
        pvalue_index,
        chr_index,
        bp_index
    )
