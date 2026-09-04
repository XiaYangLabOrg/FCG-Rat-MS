import pickle
import pandas as pd

FILE = "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SATURN/SATURN_output/saturn_results/test256_data_MS_mouse_clean_MS_FCG_rat_clean_org_saturn_seed_0_genes_to_macrogenes.pkl"

with open(FILE, "rb") as f:
    x = pickle.load(f)

rows = []

for mg in range(2000):

    tmp = []

    for gene, weights in x.items():
        tmp.append((gene, float(weights[mg])))

    tmp.sort(key=lambda z: z[1], reverse=True)

    for gene, weight in tmp[:50]:
        rows.append([mg, gene, weight])

df = pd.DataFrame(
    rows,
    columns=["macrogene", "gene", "weight"]
)

df.to_csv(
    "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SATURN/macrogene_top50_genes.csv",
    index=False
)

print(df.shape)