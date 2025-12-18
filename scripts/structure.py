import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri, conversion, default_converter
import rpy2.robjects as robjects


def seurat_to_anndata(rds_path):
    print(f"--- Loading R environment and {rds_path} ---")
    base = importr('base')
    seurat = importr('Seurat')

    robjects.globalenv['seurat_obj'] = base.readRDS(rds_path)
    seurat_obj = robjects.globalenv['seurat_obj']

    with conversion.localconverter(default_converter + pandas2ri.converter):
        # 1. Identify Assay
        active_assay = r('DefaultAssay')(seurat_obj)[0]
        print(f"Active Assay: {active_assay}")

        # 2. Extract Matrix (Sparse-aware)
        # We try 'data' first because 'integrated' assays usually lack 'counts'
        print("Extracting matrix (searching layers)...")
        r_code = f"""
        # Try data layer first, then counts
        mat <- tryCatch({{
            LayerData(seurat_obj, assay='{active_assay}', layer='data')
        }}, error = function(e) {{
            GetAssayData(seurat_obj, assay='{active_assay}', slot='data')
        }})
        if (nrow(mat) == 0 || ncol(mat) == 0) {{
             mat <- GetAssayData(seurat_obj, slot='counts')
        }}
        mat
        """
        counts_r = r(r_code)

        # Convert R sparse matrix to Scipy sparse matrix (Memory Efficient)
        # We use a manual conversion to avoid the py2rpy error
        from rpy2.robjects.vectors import Matrix
        print("Converting to Scipy sparse matrix...")
        row_names = list(r('rownames')(counts_r))
        col_names = list(r('colnames')(counts_r))

        # Use R to get the components of the sparse matrix (dgCMatrix)
        i = np.array(r('slot')(counts_r, "i"))
        p = np.array(r('slot')(counts_r, "p"))
        x = np.array(r('slot')(counts_r, "x"))
        dims = np.array(r('slot')(counts_r, "Dim"))

        X_sparse = sp.csc_matrix((x, i, p), shape=(dims[0], dims[1]))

        # 3. Extract Metadata
        metadata_py = conversion.rpy2py(r('slot')(seurat_obj, "meta.data"))

        # 4. Extract Reductions (Embeddings)
        embeddings = {}
        red_names = list(r('names')(r('slot')(seurat_obj, "reductions")))
        for name in red_names:
            try:
                # Direct array conversion to avoid dictionary conversion errors
                emb_data = np.array(r('Embeddings')(seurat_obj, reduction=name))
                embeddings[f"X_{name.lower()}"] = emb_data
                print(f" - Extracted: {name}")
            except Exception as e:
                print(f" - Skip {name}: {e}")

    # 5. Assembly
    print(f"Final Shape Check: Matrix {X_sparse.shape[1]} cells | Metadata {len(metadata_py)} cells")

    adata = ad.AnnData(
        X=X_sparse.T,
        obs=metadata_py,
        obsm=embeddings
    )
    adata.var_names = row_names

    return adata


if __name__ == "__main__":
    path = "/home/axelm@malaghan.org.nz/seurat-to-anndata/data/rds/KH_combined_2023-Jan-11.rds"
    try:
        adata = seurat_to_anndata(path)
        out = path.replace(".rds", ".h5ad")
        adata.write_h5ad(out)
        print(f"✅ Success: {out}")
    except Exception as e:
        print(f"FATAL: {e}")