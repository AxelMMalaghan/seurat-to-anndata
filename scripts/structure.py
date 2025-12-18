import os
import sys
import numpy as np
import anndata as ad
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri, conversion, default_converter

def seurat_to_anndata(rds_path):
    print(f"--- Loading R environment and {rds_path} ---")

    # 1. Initialize base R
    base = importr('base')

    # 2. FORCE R to see your user library
    # This is the "Nuclear Option" that bypasses environment variable issues
    user_lib = "/home/axelm@malaghan.org.nz/R/x86_64-pc-linux-gnu-library/4.5"
    if os.path.exists(user_lib):
        # This tells R: "Add this specific folder to your search path"
        base._libPaths(base.c(user_lib, base._libPaths()))
    else:
        print(f"Warning: User library path not found at {user_lib}")

    # 3. Now try to load Seurat
    try:
        seurat = importr('Seurat')
        print("✅ Seurat successfully loaded.")
    except Exception as e:
        print(f"❌ Failed to load Seurat even after path injection.")
        raise e

    # 4. Load the RDS file
    seurat_obj = base.readRDS(rds_path)

    with conversion.localconverter(default_converter + pandas2ri.converter):
        # Extract Assay and Counts
        active_assay = r('DefaultAssay')(seurat_obj)[0]
        print(f"Active Assay: {active_assay}")

        # Extract counts (as matrix)
        counts_r = r('as.matrix')(r('GetAssayData')(seurat_obj, slot="counts"))
        counts_py = np.array(counts_r)

        # Extract Metadata
        metadata_r = r('slot')(seurat_obj, "meta.data")
        metadata_py = conversion.rpy2py(metadata_r)

        # Extract Embeddings
        embeddings = {}
        try:
            reductions = r('names')(r('slot')(seurat_obj, "reductions"))
            for red in reductions:
                emb = r('Embeddings')(seurat_obj, reduction=red)
                embeddings[f"X_{red.lower()}"] = np.array(emb)
        except Exception:
            pass

    # 5. Build AnnData
    adata = ad.AnnData(X=counts_py.T, obs=metadata_py, obsm=embeddings)
    adata.var_names = list(r('rownames')(counts_r))
    return adata

if __name__ == "__main__":
    path_to_file = "/home/axelm@malaghan.org.nz/seurat-to-anndata/data/rds/KH_combined_2023-Jan-11.rds"
    try:
        adata = seurat_to_anndata(path_to_file)
        print(f"Successfully converted! Shape: {adata.shape}")
        output_h5ad = path_to_file.replace(".rds", ".h5ad")
        adata.write_h5ad(output_h5ad)
        print(f"Saved to: {output_h5ad}")
    except Exception as e:
        print(f"Error: {e}")