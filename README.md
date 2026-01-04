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
import logging
from core.seurat_factory import SeuratFactory

# Optional: enable logging to see progress
logging.basicConfig(level=logging.INFO)

# Load RDS and auto-detect Seurat version
extractor = SeuratFactory.load("/path/to/seurat_object.rds")

# Convert to AnnData
adata = extractor.to_anndata()

# Save
adata.write_h5ad("output.h5ad")
```

### Specify Assay and Layers

```python
# Extract normalized data as X
adata = extractor.to_anndata(assay_name="RNA")

# Extract multiple layers: first becomes X, rest go to adata.layers
adata = extractor.to_anndata(
    assay_name="RNA",
    layers=["data", "counts"]  # data -> X, counts -> adata.layers["counts"]
)

# Access layers
adata.X                    # Normalized data
adata.layers["counts"]     # Raw counts
```

### Error Handling

```python
from core.seurat_factory import SeuratFactory, SeuratLoadError

try:
    extractor = SeuratFactory.load("/path/to/file.rds")
    adata = extractor.to_anndata()
except SeuratLoadError as e:
    print(f"Failed to load Seurat object: {e}")
```

## How It Works

### 1. Load the Seurat Object

- Validates the file exists and has `.rds` extension
- Uses `rpy2` to interface with R
- Reads `.rds` files directly into R's global environment
- Validates the object is a Seurat object
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

### 3. Extract Cell Metadata (`obs`)

- Converts `seurat_obj@meta.data` to a pandas DataFrame
- Uses `pandas2ri` converter within a `localconverter` context
- Contains per-cell annotations (clusters, sample IDs, QC metrics, etc.)

### 4. Extract Gene Metadata (`var`)

- Extracts `seurat_obj@assays$<assay>@meta.features` if available
- Contains per-gene annotations (highly variable status, mean expression, etc.)
- Gracefully handles missing metadata

### 5. Extract Dimensionality Reductions

- Iterates through `seurat_obj@reductions` (UMAP, PCA, t-SNE, etc.)
- Converts each embedding matrix to numpy arrays
- Stores with `X_` prefix following AnnData conventions (`X_umap`, `X_pca`)
- Logs warnings for any failed extractions

### 6. Extract Multiple Layers

- Supports extracting multiple data layers (counts, data, scale.data)
- First layer specified becomes `adata.X`
- Additional layers stored in `adata.layers` dictionary

### 7. Assemble AnnData

- `X`: Transposed sparse matrix (Seurat is genes x cells, AnnData is cells x genes)
- `obs`: Cell metadata DataFrame
- `var`: Gene metadata DataFrame
- `obsm`: Embeddings dictionary
- `layers`: Additional expression matrices
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

## What Gets Extracted

| Component | AnnData Location | Source |
|-----------|------------------|--------|
| Expression matrix | `X` | Primary layer (default: `data`) |
| Additional layers | `layers["counts"]`, etc. | Specified in `layers` param |
| Cell metadata | `obs` | `seurat_obj@meta.data` |
| Gene metadata | `var` | `assay@meta.features` |
| Embeddings | `obsm["X_umap"]`, etc. | `seurat_obj@reductions` |
| Gene names | `var_names` | Matrix rownames |
| Cell barcodes | `obs_names` | Matrix colnames |

## Limitations

| Not Yet Supported |
|-------------------|
| Neighbor graphs (SNN/KNN) → `obsp` |
| Multiple assays simultaneously |
| Variable feature selections |

For multi-modal data (CITE-seq, multiome), run the extractor separately for each assay.

## Logging

The library uses Python's standard `logging` module. Configure it to see detailed progress:

```python
import logging

# INFO level: see major steps
logging.basicConfig(level=logging.INFO)

# DEBUG level: see all extraction details
logging.basicConfig(level=logging.DEBUG)
```

## License

MIT
