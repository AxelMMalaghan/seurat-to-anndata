import pandas as pd
import anndata as ad
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr


def load_seurat_to_anndata(path):
    # Activate the automatic conversion between R and Pandas
    pandas2ri.activate()

    # Import necessary R base package
    base = importr('base')

    print(f"Loading RDS from {path}...")
    # Read the RDS file using R's native engine
    # This avoids the "Unknown class: Seurat" error
    seurat_obj = base.readRDS(path)

    # 1. Extract Metadata (Standard R Dataframe -> Pandas)
    # Accessing the @meta.data slot
    meta_data = r('function(obj) obj@meta.data')(seurat_obj)

    # 2. Extract Count Matrix
    # We use an R snippet to pull the counts.
    # Assumes standard Seurat v3/v4/v5 structure
    print("Extracting counts...")
    counts = r('''
        function(obj) {
            library(Seurat)
            # Try to get counts from the default assay
            as.matrix(GetAssayData(obj, slot = "counts"))
        }
    ''')(seurat_obj)

    # 3. Extract PCA/UMAP Embeddings (Optional but recommended)
    # Transpose because R is Genes x Cells, AnnData is Cells x Genes
    adata = ad.AnnData(X=counts.T, obs=meta_data)

    # 4. Cleanup
    pandas2ri.deactivate()

    return adata