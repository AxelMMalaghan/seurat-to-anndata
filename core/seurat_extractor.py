import logging
from abc import ABC, abstractmethod
from typing import List, Optional

import anndata as ad
import numpy as np
from rpy2.robjects import r, pandas2ri, conversion, default_converter

logger = logging.getLogger(__name__)


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

    def to_anndata(
        self,
        assay_name: Optional[str] = None,
        layers: Optional[List[str]] = None
    ) -> ad.AnnData:
        """
        Orchestrates the decomposition of the Seurat object.

        Args:
            assay_name: Assay to extract (None = default assay)
            layers: List of layers to extract (e.g., ['counts', 'data']).
                    First layer becomes X, rest go to adata.layers.
                    Defaults to ['data'] if None.

        Returns:
            Fully assembled AnnData object
        """
        if assay_name is None:
            assay_name = r(f'DefaultAssay({self.r_var_name})')[0]
            logger.info(f"Using default assay: {assay_name}")

        if layers is None:
            layers = ["data"]

        # Extract first layer as X
        primary_layer = layers[0]
        logger.info(f"Extracting primary layer '{primary_layer}' from assay '{assay_name}'")
        X_sparse, var_names, obs_names = self._extract_assay(assay_name, primary_layer)

        obs = self._extract_obs()
        var = self._extract_var(assay_name)

        adata = ad.AnnData(
            X=X_sparse.T,  # Transpose: Seurat is genes x cells, AnnData is cells x genes
            obs=obs,
            var=var
        )
        adata.var_names = var_names
        adata.obs_names = obs_names

        # Extract additional layers
        for layer_name in layers[1:]:
            self._add_layer(adata, assay_name, layer_name)

        self._add_embeddings(adata)

        logger.info(f"AnnData created: {adata.n_obs} cells x {adata.n_vars} genes")
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

    @abstractmethod
    def _add_layer(self, adata: ad.AnnData, assay_name: str, layer_name: str) -> None:
        """
        Extract an additional layer and add to adata.layers.

        Must be implemented in subclasses.
        """
        pass

    def _extract_obs(self):
        """
        Extract cell metadata from seurat_obj@meta.data.

        Returns:
            pandas DataFrame with cell metadata
        """
        logger.debug("Extracting cell metadata (obs)")
        with conversion.localconverter(default_converter + pandas2ri.converter):
            obs = conversion.rpy2py(r(f'{self.r_var_name}@meta.data'))
        logger.debug(f"Extracted {len(obs)} cells with {len(obs.columns)} metadata columns")
        return obs

    def _extract_var(self, assay_name: str):
        """
        Extract gene/feature metadata from the assay's meta.features.

        Args:
            assay_name: Name of the assay to extract metadata from

        Returns:
            pandas DataFrame with gene metadata (may be empty)
        """
        logger.debug(f"Extracting gene metadata (var) from assay '{assay_name}'")
        try:
            # Check if meta.features exists and has data
            has_meta = r(f'!is.null({self.r_var_name}@assays${assay_name}@meta.features) && nrow({self.r_var_name}@assays${assay_name}@meta.features) > 0')[0]
            if has_meta:
                with conversion.localconverter(default_converter + pandas2ri.converter):
                    var = conversion.rpy2py(r(f'{self.r_var_name}@assays${assay_name}@meta.features'))
                logger.debug(f"Extracted {len(var.columns)} gene metadata columns")
                return var
            else:
                logger.debug("No gene metadata found, returning empty DataFrame")
                return None
        except Exception as e:
            logger.warning(f"Could not extract gene metadata: {e}")
            return None

    def _add_embeddings(self, adata: ad.AnnData) -> None:
        """
        Extract reductions and add to adata.obsm in place.
        """
        red_names_r = r(f'names({self.r_var_name}@reductions)')
        if red_names_r == r('NULL'):
            logger.debug("No reductions found in Seurat object")
            return

        red_names = list(red_names_r)
        logger.debug(f"Found {len(red_names)} reductions: {red_names}")

        for name in red_names:
            try:
                emb = np.array(r(f'as.matrix(Embeddings({self.r_var_name}, reduction="{name}"))'))
                key = f"X_{name.lower()}"
                adata.obsm[key] = emb
                logger.debug(f"Added embedding '{key}' with shape {emb.shape}")
            except Exception as e:
                logger.warning(f"Failed to extract reduction '{name}': {e}")
