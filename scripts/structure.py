import os
import numpy as np
import anndata as ad
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri, conversion, default_converter


def seurat_to_anndata(rds_path):
    print(f"--- Loading R environment and {rds_path} ---")

    # No manual path injection here!
    base = importr('base')

    try:
        seurat = importr('Seurat')
        print("✅ Seurat successfully loaded from Conda environment.")
    except Exception as e:
        print(f"❌ Seurat not found in Conda env. Run 'conda install r-seurat'")
        raise e

    seurat_obj = base.readRDS(rds_path)

    with conversion.localconverter(default_converter + pandas2ri.converter):
        active_assay = r('DefaultAssay')(seurat_obj)[0]
        # Use GetAssayData to extract the counts
        counts_r = r('as.matrix')(r('GetAssayData')(seurat_obj, slot="counts"))
        counts_py = np.array(counts_r)
        metadata_py = conversion.rpy2py(r('slot')(seurat_obj, "meta.data"))

        embeddings = {}
        try:
            reductions = r('names')(r('slot')(seurat_obj, "reductions"))
            for red in reductions:
                emb = r('Embeddings')(seurat_obj, reduction=red)
                embeddings[f"X_{red.lower()}"] = np.array(emb)
        except Exception:
            pass

    adata = ad.AnnData(X=counts_py.T, obs=metadata_py, obsm=embeddings)
    adata.var_names = list(r('rownames')(counts_r))
    return adata


if __name__ == "__main__":
    path_to_file = "/home/axelm@malaghan.org.nz/seurat-to-anndata/data/rds/KH_combined_2023-Jan-11.rds"
    adata = seurat_to_anndata(path_to_file)
    output_h5ad = path_to_file.replace(".rds", ".h5ad")
    adata.write_h5ad(output_h5ad)
    print(f"Done! Saved to {output_h5ad}")