# seurat-to-anndata

Convert Seurat objects (.rds) to AnnData (.h5ad) format for use with Scanpy and the Python single-cell ecosystem.

## Installation

### Requirements

- Python 3.8+
- R 4.0+ with Seurat installed
- rpy2

```bash
pip install anndata numpy scipy rpy2 pandas
```

Ensure R and Seurat are available:

```r
install.packages("Seurat")
```

## Usage

### Quick Start

```python
from core.seurat_factory import SeuratFactory

# Load RDS and auto-detect Seurat version
extractor = SeuratFactory.load("/path/to/seurat_object.rds")

# Convert to AnnData
adata = extractor.to_anndata()

# Save
adata.write_h5ad("output.h5ad")
```

### Specify Assay and Layer

```python
# Extract raw counts from RNA assay
adata = extractor.to_anndata(assay_name="RNA", layer_name="counts")

# Extract normalized data from integrated assay
adata = extractor.to_anndata(assay_name="integrated", layer_name="data")
```

## How It Works

### 1. Load the Seurat Object

- Uses `rpy2` to interface with R
- Reads `.rds` files directly into R's global environment
- Auto-detects Seurat version (v4 vs v5) based on assay class

### 2. Extract Expression Matrix

The primary challenge is transferring R's sparse matrix to Python. The converter:

- Uses `LayerData()` for Seurat v5, falls back to `GetAssayData()` for v4
- Extracts raw dgCMatrix components individually to avoid rpy2 conversion errors:
  - `i`: row indices
  - `p`: column pointers
  - `x`: non-zero values
  - `Dim`: matrix dimensions
- Reconstructs a `scipy.sparse.csc_matrix` from these primitives

### 3. Extract Cell Metadata

- Converts `seurat_obj@meta.data` to a pandas DataFrame
- Uses `pandas2ri` converter within a `localconverter` context
- Contains per-cell annotations (clusters, sample IDs, QC metrics, etc.)

### 4. Extract Dimensionality Reductions

- Iterates through `seurat_obj@reductions` (UMAP, PCA, t-SNE, etc.)
- Converts each embedding matrix to numpy arrays
- Stores with `X_` prefix following AnnData conventions (`X_umap`, `X_pca`)

### 5. Assemble AnnData

- `X`: Transposed sparse matrix (Seurat is genes x cells, AnnData is cells x genes)
- `obs`: Cell metadata DataFrame
- `obsm`: Embeddings dictionary
- `var_names`: Gene names
- `obs_names`: Cell barcodes

## Architecture

```
seurat-to-anndata/
├── core/
│   ├── seurat_extractor.py   # Abstract base class
│   ├── seurat_factory.py     # Loads RDS, returns version-specific extractor
│   └── __init__.py
├── seurat_v4_extractor.py    # Seurat v4 (GetAssayData API)
├── seurat_v5_extractor.py    # Seurat v5 (LayerData API)
└── scripts/
    └── working.py            # Standalone conversion script
```

## Limitations

| Extracted | Not Extracted |
|-----------|---------------|
| Expression matrix (one layer) | Multiple layers simultaneously |
| Cell metadata (`obs`) | Gene metadata (`var`) |
| Embeddings (UMAP, PCA, etc.) | Neighbor graphs (SNN/KNN) |
| Active assay | Additional assays (ADT, ATAC) |

For multi-modal data or full fidelity conversion, additional extraction logic would be needed.

## License

MIT
