import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri, conversion, default_converter
import rpy2.robjects as robjects


def seurat_to_anndata(rds_path):
    print(f"--- Loading R environment and {rds_path} ---")

    # 1. Initialize R libraries
    base = importr('base')
    try:
        seurat = importr('Seurat')
        print("✅ Seurat successfully loaded.")
    except Exception as e:
        print(f"❌ Seurat not found in Conda env.")
        raise e

    # 2. Load the RDS file into R
    # We assign it to a named variable in R to make indexing easier
    robjects.globalenv['seurat_obj'] = base.readRDS(rds_path)
    seurat_obj = robjects.globalenv['seurat_obj']

    with conversion.localconverter(default_converter + pandas2ri.converter):
        # 3. Identify the Active Assay
        active_assay = r('DefaultAssay')(seurat_obj)[0]
        print(f"Active Assay detected: {active_assay}")

        # 4. Extract Matrix (v4 and v5 compatible)
        # We try LayerData (v5) first, then GetAssayData (v4)
        print("Extracting count matrix...")
        r_code = f"""
        tryCatch({{
            as.matrix(LayerData(seurat_obj, assay='{active_assay}', layer='counts'))
        }}, error = function(e) {{
            as.matrix(GetAssayData(seurat_obj, assay='{active_assay}', slot='counts'))
        }})
        """
        counts_mat = r(r_code)
        counts_py = np.array(counts_mat)

        # 5. Extract Metadata (obs)
        print("Extracting cell metadata (obs)...")
        metadata_r = r('slot')(seurat_obj, "meta.data")
        metadata_py = conversion.rpy2py(metadata_r)

        # 6. Extract Dimensional Reductions (obsm)
        print("Extracting embeddings (obsm)...")
        embeddings = {}
        try:
            reductions = r('names')(r('slot')(seurat_obj, "reductions"))
            for red in reductions:
                emb = r('Embeddings')(seurat_obj, reduction=red)
                # Map names to scanpy standard (e.g., X_umap)
                key = f"X_{red.lower()}"
                embeddings[key] = np.array(emb)
                print(f" - Found embedding: {red}")
        except Exception as e:
            print(f" - Note: Could not extract reductions: {e}")

    # 7. Final Validation & Assembly
    print(f"Final Check: Matrix has {counts_py.shape[1]} cells. Metadata has {len(metadata_py)} cells.")

    if counts_py.shape[1] != len(metadata_py):
        print("⚠️ Warning: Cell count mismatch! Attempting to subset metadata to match matrix...")
        # Subset metadata based on the column names of the R matrix
        cell_names = list(r('colnames')(counts_mat))
        metadata_py = metadata_py.loc[cell_names]

    # Transpose matrix for AnnData (Cells x Genes)
    adata = ad.AnnData(
        X=counts_py.T,
        obs=metadata_py,
        obsm=embeddings
    )

    # Set Gene Names (var)
    adata.var_names = list(r('rownames')(counts_mat))

    return adata


if __name__ == "__main__":
    # Path configuration
    input_rds = "/home/axelm@malaghan.org.nz/seurat-to-anndata/data/rds/KH_combined_2023-Jan-11.rds"
    output_h5ad = input_rds.replace(".rds", ".h5ad")

    try:
        adata_obj = seurat_to_anndata(input_rds)

        print("\n--- Success! ---")
        print(adata_obj)

        # Writing to disk
        print(f"Saving AnnData to {output_h5ad}...")
        adata_obj.write_h5ad(output_h5ad)
        print("File saved successfully.")

    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        sys.exit(1)