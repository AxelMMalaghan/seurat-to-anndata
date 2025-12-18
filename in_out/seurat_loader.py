import os
import subprocess
import pandas as pd
import anndata as ad
from scipy.sparse import csc_matrix


class SeuratLoader:
    def __init__(self):
        self._setup_r_environment()
        import rpy2.robjects as robjects
        from rpy2.robjects import pandas2ri
        self.robjects = robjects
        self.r = robjects.r
        pandas2ri.activate()

    def _setup_r_environment(self):
        try:
            r_ver_out = subprocess.check_output(["R", "--version"]).decode()
            r_ver = r_ver_out.splitlines()[0].split()[2]
            major_minor = ".".join(r_ver.split(".")[:2])
            lib_path = os.path.expanduser(f"~/R/x86_64-pc-linux-gnu-library/{major_minor}")
            os.makedirs(lib_path, exist_ok=True)
            os.environ['R_LIBS_USER'] = lib_path
        except Exception as e:
            raise RuntimeError(f"Could not configure R environment: {e}")

    def inspect_rds(self, path):
        """Quickly inspects the RDS without full conversion."""
        print(f"--- Inspecting: {os.path.basename(path)} ---")
        self.r(f'seurat_obj <- readRDS("{path}")')

        info = {
            "version": str(self.r('as.character(seurat_obj@version)')[0]),
            "assays": list(self.r('names(seurat_obj@assays)')),
            "active_assay": str(self.r('Seurat::DefaultAssay(seurat_obj)')[0]),
            "reductions": list(self.r('names(seurat_obj@reductions)')),
            "cell_count": int(self.r('ncol(seurat_obj)')[0]),
            "gene_count": int(self.r('nrow(seurat_obj)')[0])
        }

        print(f"Found Seurat v{info['version']} with {info['cell_count']} cells.")
        return info

    def load_data(self, path):
        if not os.path.exists(path):
            raise FileNotFoundError(f"RDS file not found: {path}")

        try:
            # 1. Inspection & Setup
            metadata = self.inspect_rds(path)
            is_v5 = metadata['version'].startswith('5')
            default_assay = metadata['active_assay']

            # 2. Extract Metadata (obs)
            obs = self.r('seurat_obj@meta.data')

            # 3. Extract Gene Metadata (var) & Variable Features
            genes = list(self.r('rownames(seurat_obj)'))
            var = pd.DataFrame(index=genes)

            # Sync Variable Features
            hvgs = list(self.r(f'Seurat::VariableFeatures(seurat_obj, assay = "{default_assay}")'))
            var['highly_variable'] = var.index.isin(hvgs)

            # 4. Extract Sparse Matrix (X)
            r_getter = "SeuratObject::LayerData" if is_v5 else "Seurat::GetAssayData"
            X = self._get_sparse_matrix(f"{r_getter}(seurat_obj, assay = '{default_assay}', slot = 'counts')")

            # 5. Build AnnData
            adata = ad.AnnData(X=X.T, obs=obs, var=var)
            adata = self._add_embeddings(adata)

            return adata

        except Exception as e:
            raise IOError(f"Seurat RDS file could not be loaded: {e}")
        finally:
            self.r('rm(seurat_obj); gc()')

    def _get_sparse_matrix(self, r_expression):
        sparse_parts = self.r(f'''
            library(Matrix)
            m <- {r_expression}
            list(data=m@x, i=m@i, p=m@p, dim=dim(m))
        ''')
        return csc_matrix(
            (sparse_parts.rx2('data'), sparse_parts.rx2('i'), sparse_parts.rx2('p')),
            shape=tuple(sparse_parts.rx2('dim'))
        )

    def _add_embeddings(self, adata):
        reductions = list(self.r('names(seurat_obj@reductions)'))
        for name in reductions:
            try:
                coords = self.r(f'as.data.frame(Seurat::Embeddings(seurat_obj, reduction = "{name}"))')
                adata.obsm[f'X_{name.lower()}'] = coords.values
            except:
                continue
        return adata