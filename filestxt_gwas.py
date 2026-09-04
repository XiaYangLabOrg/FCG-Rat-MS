# %%
import os

# ---- USER SETTINGS ----
gwas_dir = '/u/scratch/v/vturnbil/GWAS'
output_file = '/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/Mergeomics/MDF/newfiles.txt'
extensions = ('.tsv', '.txt')

# Column name variants (lowercase)
CHR_NAMES = {'chr', 'chromosome'}
BP_NAMES  = {'bp', 'pos', 'position', 'base_pair', 'base_pair_location'}
P_NAMES   = {'p', 'pval', 'pvalue', 'p_value', 'p-value'}

def detect_columns(header):
    chr_idx = bp_idx = p_idx = None

    for i, col in enumerate(header):
        col_norm = col.lower()

        if chr_idx is None and col_norm in CHR_NAMES:
            chr_idx = i
        elif bp_idx is None and col_norm in BP_NAMES:
            bp_idx = i
        elif p_idx is None and col_norm in P_NAMES:
            p_idx = i

    return chr_idx, bp_idx, p_idx


written = []
skipped = []

print("\nGWAS column detection summary")
print("=" * 60)

# ---- Scan all GWAS files ----
for fname in sorted(os.listdir(gwas_dir)):
    if not fname.endswith(extensions):
        continue

    path = os.path.join(gwas_dir, fname)
    base = os.path.splitext(fname)[0]

    # Read first non-empty line (header)
    try:
        with open(path) as f:
            for line in f:
                if line.strip():
                    header = line.strip().split()
                    break
            else:
                skipped.append((base, 'EMPTY FILE'))
                continue
    except Exception as e:
        skipped.append((base, f'UNREADABLE FILE: {e}'))
        continue

    chr_idx, bp_idx, p_idx = detect_columns(header)

    # ---- Missing column checks ----
    if chr_idx is None:
        skipped.append((base, 'MISSING CHROMOSOME COLUMN'))
        continue

    if bp_idx is None:
        skipped.append((base, 'MISSING BASEPAIR COLUMN'))
        continue

    if p_idx is None:
        skipped.append((base, 'MISSING P-VALUE COLUMN'))
        continue

    # ---- Record entry ----
    written.append(f'{base}({chr_idx},{bp_idx},{p_idx})')

    # ---- Print detected columns for verification ----
    print(f'\n{base}')
    print('-' * len(base))
    print(f'Chromosome: index {chr_idx} → "{header[chr_idx]}"')
    print(f'Basepair:   index {bp_idx} → "{header[bp_idx]}"')
    print(f'P-value:    index {p_idx} → "{header[p_idx]}"')


# ---- Write newfiles.txt ----
with open(output_file, 'w') as out:
    for line in written:
        out.write(line + '\n')


# ---- Final report ----
print('\nnewfiles.txt generation report')
print('=' * 60)

print(f'Written entries: {len(written)}')
print(f'Skipped entries: {len(skipped)}')

if skipped:
    print('\nSkipped files (with reasons):')
    for name, reason in skipped:
        print(f'  {name}: {reason}')
else:
    print('\nNo files were skipped.')
