import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csc_matrix
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.packages import importr, isinstalled


class RobustSeuratLoader:
    def __init__(self):
        # 1. Ensure R environment is ready
        self._check_r_dependencies()
        pandas2ri.activate()

    def _check_r_dependencies(self):
        required_packages = ['Seurat', 'Matrix', 'base']
        for pkg in required_packages:
            if not isinstalled(pkg):
                raise ImportError(f"R package '{pkg}' is not installed. "
                                  "Please run 'install.packages(\"{pkg}\")' in R.")

    def load_data(self, path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"RDS file not found at: {path}")

        print(f"--- Loading RDS: {os.path.basename(path)} ---")

        # Load the object into the R global environment
        robjects.r(f'seurat_obj <- readRDS("{path}")')

        # Determine Seurat version to adjust extraction logic
        s_version = robjects.r('as.character(seurat_obj@version)')[0]
        print(f"Detected Seurat version: {s_version}")

        try:
            # 1. EXTRACT METADATA
            obs = self._extract_metadata()

            # 2. EXTRACT COUNTS (Handling Sparsity)
            X = self._extract_sparse_counts()

            # 3. EXTRACT VARIABLE FEATURES (If they exist)
            var = self._extract_variable_genes()

            # Create AnnData
            adata = ad.AnnData(X=X, obs=obs, var=var)

            # 4. EXTRACT EMBEDDINGS (UMAP, PCA)
            adata = self._extract_embeddings(adata)

            print("Successfully converted to AnnData.")
            return adata

        except Exception as e:
            raise RuntimeError(f"Failed to parse Seurat object: {str(e)}")
        finally:
            # Clean up R memory
            robjects.r('rm(seurat_obj); gc()')

    def _extract_metadata(self):
        print("Extracting metadata...")
        return robjects.r('seurat_obj@meta.data')

    def _extract_sparse_counts(self):
        """
        Extracts counts as a Scipy Sparse Matrix to prevent MemoryError
        """
        print("Extracting sparse counts (this may take a moment)...")
        # We pull the components of the sparse matrix (i, p, x) from R
        # This is much faster and memory-safe than converting to a dense array
        r_script = """
        library(Matrix)
        counts <- Seurat::GetAssayData(seurat_obj, slot = "counts")
        list(
            data = counts@x,
            indices = counts@i,
            indptr = counts@p,
            dim = dim(counts),
            row_names = rownames(counts),
            col_names = colnames(counts)
        )
        """
        sparse_data = robjects.r(r_script)

        # Reconstruct csc_matrix in Python
        # R is 0-indexed for 'i' and 'p' in newer Matrix versions,
        # but rpy2 handles the conversion gracefully.
        mtx = csc_matrix(
            (sparse_data.rx2('data'), sparse_data.rx2('indices'), sparse_data.rx2('indptr')),
            shape=tuple(sparse_data.rx2('dim'))
        )
        # AnnData expects (cells x genes), Seurat is (genes x cells)
        return mtx.T

    def _extract_variable_genes(self):
        gene_names = robjects.r('rownames(seurat_obj)')
        var_df = pd.DataFrame(index=list(gene_names))

        # Mark highly variable genes if present
        hvg = robjects.r('VariableFeatures(seurat_obj)')
        if hvg is not robjects.NULL:
            var_df['highly_variable'] = var_df.index.isin(list(hvg))
        return var_df

    def _extract_embeddings(self, adata):
        embeddings = robjects.r('names(seurat_obj@reductions)')
        if embeddings is robjects.NULL:
            return adata

        for emb in list(embeddings):
            print(f"Extracting embedding: {emb}")
            coord_script = f'as.data.frame(Embeddings(seurat_obj, reduction = "{emb}"))'
            coords = robjects.r(coord_script)
            adata.obsm[f'X_{emb.lower()}'] = coords.values
        return adata