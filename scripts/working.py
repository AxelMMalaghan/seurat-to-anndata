import numpy as np
import anndata as ad
import scipy.sparse as sp
from rpy2.robjects import r, pandas2ri, conversion, default_converter


"""
How this works:

1) Load the Seurat Object

    - Imports R packages via rpy2
    - Reads .rds files directly into R's global env
    - Identifies the active assay (e.g "RNA", "Integrated")
    
2) Extract Expression Matrix

    The primary challenge of this is getting the sparse matrix from R to Python. 
    
    The script:
    
        - Attempts to get data via LayerData() (Seurat v5), before falling back to GetAssayData() (v4)
        - Extracts the raw dgCMatrix components individually
        
            - i: row indices
            - p: col pointers
            - x: non-zero values
            - Dim: matrix dimensions 
            
        - Reconstructs a scipy.sparse.csc_matrix from these primitives

    The manual extraction avoids rpy2's OrderedDict conversion errors that occur with some R objects.
    
3) Extract Cell Metadata

    - Uses pandas2ri.converter within a localconverter context
    - Directly converts seurat_obj@meta.data to a pandas DataFrame
    - This contains the per-cell annoations
    
4) Extract Dimensionality Reductions

    - Iterates through all reductions in seurat_obj@reductions (UMAP, PCA, TSNE, etc)
    - Converts each embedding matrix to numpy arrays
    - Stores them with X_ prefix (AnnData Conversion)
    
5) Assemble AnnData

    - Creates AnnData with:
    
        - X: transposed sparse matrix (R is wrong way round (m x n))
        - obs: cell metadata DataFrame
        - obsm: embeddings dictionary
        
    - Sets var_names (genes) and obs_names (cells)
    
Output

Hopefully the resulting file contains the underlying data and maybe even in the correct way!
"""





def seurat_to_anndata(rds_path):

    print(f"--- Loading R environment and {rds_path} ---")

    # Load into R global environment
    r(f'seurat_obj <- readRDS("{rds_path}")')

    # 1. Identify Assay
    active_assay = r('DefaultAssay(seurat_obj)')[0]
    print(f"Active Assay: {active_assay}")

    # 2. Extract Sparse Matrix Components manually as raw NumPy arrays
    print("Extracting sparse matrix components...")
    r('''
        # Prioritize integrated data layer, fallback to RNA counts
        mat <- tryCatch({
            LayerData(seurat_obj, assay=''' + f"'{active_assay}'" + ''', layer='data')
        }, error = function(e) {
            GetAssayData(seurat_obj, slot='data')
        })
    ''')

    # Extract primitives to bypass OrdDict errors
    i = np.array(r('as.integer(mat@i)'))
    p = np.array(r('as.integer(mat@p)'))
    x = np.array(r('as.numeric(mat@x)'))
    dims = np.array(r('as.integer(mat@Dim)'))
    row_names = list(r('rownames(mat)'))
    col_names = list(r('colnames(mat)'))

    X_sparse = sp.csc_matrix((x, i, p), shape=(dims[0], dims[1]))
    print(f"Matrix extracted: {X_sparse.shape[0]} genes, {X_sparse.shape[1]} cells")

    # 3. Extract Metadata (Converting to CSV in R memory to avoid OrdDict)
    # This is a bulletproof way to get metadata into Pandas
    print("Extracting metadata via CSV bridge...")
    metadata_csv = r(
        'write.table(seurat_obj@meta.data, sep=",", col.names=NA, quote=FALSE, file=rawConnection(raw(0), "w"))')
    # Use a simpler approach: direct conversion within the localconverter
    with conversion.localconverter(default_converter + pandas2ri.converter):
        metadata_py = conversion.rpy2py(r('seurat_obj@meta.data'))

    # 4. Extract Reductions
    embeddings = {}
    red_names = list(r('names(seurat_obj@reductions)'))
    for name in red_names:
        try:
            emb_data = np.array(r(f'as.matrix(Embeddings(seurat_obj, reduction="{name}"))'))
            embeddings[f"X_{name.lower()}"] = emb_data
            print(f" - Extracted: {name}")
        except Exception:
            pass

    # 5. Assembly
    adata = ad.AnnData(
        X=X_sparse.T,
        obs=metadata_py,
        obsm=embeddings
    )
    adata.var_names = row_names
    adata.obs_names = col_names

    return adata


if __name__ == "__main__":
    path = "/home/axelm@malaghan.org.nz/seurat-to-anndata/data/rds/KH_combined_2023-Jan-11.rds"
    try:
        adata = seurat_to_anndata(path)
        out = path.replace(".rds", ".h5ad")
        adata.write_h5ad(out)
        print(f"✅ Success! AnnData saved to: {out}")
    except Exception as e:
        print(f"FATAL: {e}")