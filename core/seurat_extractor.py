from abc import ABC, abstractmethod

import anndata as ad
import numpy as np
from rpy2.robjects import r, pandas2ri, conversion, default_converter


class SeuratExtractor(ABC):
    """
    Base interface for all Seurat extractors (based on version).
    Subclass and implement _extract_assay for version-specific logic.
    """

    def __init__(self, r_var_name: str = "seurat_obj"):
        """
        Args:
            r_var_name: Name of the Seurat object in R's global environment
        """
        self.r_var_name = r_var_name

    def to_anndata(self, assay_name: str = None, layer_name: str = "data") -> ad.AnnData:
        """
        Orchestrates the decomposition of the Seurat object by delegating to:
        - _extract_obs
        - _extract_assay
        - _add_embeddings

        Args:
            assay_name: Assay to extract (None = default assay)
            layer_name: Layer to extract ('data', 'counts', 'scale.data')

        Returns:
            Fully assembled AnnData object
        """
        if assay_name is None:
            assay_name = r(f'DefaultAssay({self.r_var_name})')[0]

        X_sparse, var_names, obs_names = self._extract_assay(assay_name, layer_name)
        obs = self._extract_obs()

        adata = ad.AnnData(
            X=X_sparse.T,  # Transpose: Seurat is genes x cells, AnnData is cells x genes
            obs=obs
        )
        adata.var_names = var_names
        adata.obs_names = obs_names

        self._add_embeddings(adata)

        return adata

    @abstractmethod
    def _extract_assay(self, assay_name: str, layer_name: str) -> tuple:
        """
        Extract sparse matrix from the specified assay/layer.

        Must be implemented in subclasses - Seurat v4 vs v5 differ here.

        Returns:
            Tuple of (X_sparse, var_names, obs_names)
        """
        pass

    def _extract_obs(self):
        """
        Extract cell metadata from seurat_obj@meta.data.

        Returns:
            pandas DataFrame with cell metadata
        """
        with conversion.localconverter(default_converter + pandas2ri.converter):
            obs = conversion.rpy2py(r(f'{self.r_var_name}@meta.data'))
        return obs

    def _add_embeddings(self, adata: ad.AnnData) -> None:
        """
        Extract reductions and add to adata.obsm in place.
        """
        red_names = list(r(f'names({self.r_var_name}@reductions)'))

        for name in red_names:
            try:
                emb = np.array(r(f'as.matrix(Embeddings({self.r_var_name}, reduction="{name}"))'))
                adata.obsm[f"X_{name.lower()}"] = emb
            except Exception:
                pass  # Skip failed extractions silently
