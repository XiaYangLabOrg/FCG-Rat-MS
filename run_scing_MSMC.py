import os
import sys
from pathlib import Path
import warnings
import scanpy as sc
from scing import supercells, build, merge
import numpy as np


def main():

    # -------------------------
    # Get input file from command line
    # -------------------------
    if len(sys.argv) < 2:
        raise ValueError("Usage: python run_scing_MSMC.py <input.h5ad>")

    input_path = Path(sys.argv[1]).resolve()

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    dataset = input_path.stem   # e.g. MSMC_Astrocytes_XXO

    print(f"\n=== Running SCING on {dataset} ===")
    print(f"Input file: {input_path}\n")

    save = True
    all_edges = []

    # -------------------------
    # Load AnnData
    # -------------------------
    adata = sc.read_h5ad(input_path)

    # -------------------------
    # Output directory (auto per dataset)
    # -------------------------
    base_out = Path(
        "/u/project/xyang123/xyang123-NOBACKUP/turnbill/GSU_FCG/Restart/MS/SCING/SCING_output"
    )

    outdir = base_out / dataset
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Output directory: {outdir}\n")

    # -------------------------
    # Build GRNs
    # -------------------------
    for i in range(10):
        print(f'Building GRN {i}')

        grn = build.grnBuilder(
            adata=adata,
            ngenes=-1,
            nneighbors=100,
            npcs=10,
            subsample_perc=0.7,
            prefix=f'net.{i}',
            outdir=str(outdir),
            ncore=8,
            mem_per_core=int(3e9),
            verbose=True,
        )

        grn.subsample_cells()
        grn.filter_genes()
        grn.filter_gene_connectivities()
        grn.build_grn()

        if save:
            grn.save_edges()

        all_edges.append(grn.edges)

    # -------------------------
    # Merge networks
    # -------------------------
    print('Merging networks...')

    merger = merge.NetworkMerger(
        adata=adata,
        networks=all_edges,
        cycles=None,
        minimum_edge_appearance_threshold=0.2,
        prefix='final',
        outdir=str(outdir),
        ncore=8,
        mem_per_core=int(3e9),
        verbose=True
    )

    merger.preprocess_network_files()
    merger.remove_reversed_edges()
    merger.remove_cycles()
    merger.get_triads()
    merger.remove_redundant_edges()
    merger.save_network()

    # -------------------------
    # Sort edges by importance
    # -------------------------
    merger.edge_df = merger.edge_df.sort_values(
        by='importance',
        ascending=False
    )

    print('\nPipeline finished successfully.')


if __name__ == "__main__":
    main()
