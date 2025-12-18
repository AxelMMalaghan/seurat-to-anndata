import os
import sys
import numpy as np
import anndata as ad
# Set the R library path so rpy2 can find Seurat
R_LIB_PATH = "/home/axelm@malaghan.org.nz/R/x86_64-pc-linux-gnu-library/4.5"
os.environ['R_LIBS_USER'] = R_LIB_PATH

from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri, conversion, default_converter
import rpy2.robjects as robjects


# Define the converter globally or within the function
def seurat_to_anndata(rds_path):
    print(f"--- Loading R environment and {rds_path} ---")

    base = importr('base')
    seurat = importr('Seurat')

    # 1. Load the RDS file
    seurat_obj = base.readRDS(rds_path)

    # Use the context manager to handle conversion safely
    with conversion.localconverter(default_converter + pandas2ri.converter):

        # 2. Extract Assay and Counts
        active_assay = r('DefaultAssay')(seurat_obj)[0]
        print(f"Active Assay: {active_assay}")

        # Get counts as a matrix
        # Note: If memory is an issue, we can try to keep it sparse later
        counts_r = r('as.matrix')(r('GetAssayData')(seurat_obj, slot="counts"))
        counts_py = np.array(counts_r)

        # 3. Extract Metadata (obs)
        metadata_r = r('slot')(seurat_obj, "meta.data")
        metadata_py = conversion.rpy2py(metadata_r)

        # 4. Extract Embeddings (obsm)
        embeddings = {}
        try:
            reductions = r('names')(r('slot')(seurat_obj, "reductions"))
            for red in reductions:
                emb = r('Embeddings')(seurat_obj, reduction=red)
                embeddings[f"X_{red.lower()}"] = np.array(emb)
                print(f"Extracted embedding: {red}")
        except Exception:
            print("No embeddings found.")

    # 5. Build AnnData (Outside the converter context for better performance)
    adata = ad.AnnData(
        X=counts_py.T,
        obs=metadata_py,
        obsm=embeddings
    )

    # Set gene names (var names)
    adata.var_names = list(r('rownames')(counts_r))

    return adata


if __name__ == "__main__":
    path_to_file = "/home/axelm@malaghan.org.nz/seurat-to-anndata/data/rds/KH_combined_2023-Jan-11.rds"
    try:
        adata = seurat_to_anndata(path_to_file)
        print(f"\nSuccessfully converted! Shape: {adata.shape}")

        output_h5ad = path_to_file.replace(".rds", ".h5ad")
        adata.write_h5ad(output_h5ad)
    except Exception as e:
        print(f"Error: {e}")

