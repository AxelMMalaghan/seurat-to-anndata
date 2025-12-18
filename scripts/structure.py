import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr

# 1. Enable R to Pandas conversion
pandas2ri.activate()


def seurat_to_anndata(rds_path):
    print(f"--- Loading R environment and {rds_path} ---")

    # Import R libraries
    base = importr('base')
    seurat = importr('Seurat')

    # Load the RDS file into R memory
    seurat_obj = base.readRDS(rds_path)

    # 2. Extract the Active Assay (e.g., RNA or SCT)
    # We use R's accessor methods to ensure we get the right slot
    active_assay = r('DefaultAssay')(seurat_obj)[0]
    print(f"Active Assay detected: {active_assay}")

    # 3. Extract the Count Matrix
    # This pulls the 'counts' slot from the active assay
    # We use GetAssayData to handle sparse or dense formats correctly
    counts = r('as.matrix')(r('GetAssayData')(seurat_obj, slot="counts"))

    # 4. Extract Metadata (obs)
    # Seurat metadata is stored in the @meta.data slot
    metadata = r('as.data.frame')(r('slot')(seurat_obj, "meta.data"))

    # 5. Extract Feature Metadata (var)
    # This gets gene-level information (like gene names/IDs)
    features = r('as.data.frame')(r('slot')(r('slot')(seurat_obj, "assays").rx2(active_assay), "meta.features"))

    # 6. Extract Embeddings (obsm) - Optional but useful for browsers
    # Tries to pull UMAP or TSNE if they exist
    embeddings = {}
    try:
        reductions = r('names')(r('slot')(seurat_obj, "reductions"))
        for red in reductions:
            emb = r('Embeddings')(seurat_obj, reduction=red)
            embeddings[f"X_{red.lower()}"] = np.array(emb)
            print(f"Extracted embedding: {red}")
    except Exception as e:
        print(f"No embeddings found: {e}")

    # 7. Assemble AnnData
    # Transpose counts (Seurat is Features x Cells, AnnData is Cells x Features)
    adata = ad.AnnData(
        X=counts.T,
        obs=metadata,
        var=features if not features.empty else pd.DataFrame(index=r('rownames')(counts)),
        obsm=embeddings
    )

    return adata


# --- Execution ---
if __name__ == "__main__":
    path_to_file = "/home/axelm@malaghan.org.nz/seurat-to-anndata/data/rds/KH_combined_2023-Jan-11.rds"

    try:
        adata = seurat_to_anndata(path_to_file)
        print("\n--- Success! ---")
        print(adata)

        # Save for your browser project
        output_h5ad = path_to_file.replace(".rds", ".h5ad")
        adata.write_h5ad(output_h5ad)
        print(f"Saved to: {output_h5ad}")

    except Exception as e:
        print(f"Conversion failed: {e}")