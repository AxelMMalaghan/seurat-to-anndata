import logging

import anndata as ad
import numpy as np
import scipy.sparse as sp
from rpy2.robjects import r

from core.seurat_extractor import SeuratExtractor

logger = logging.getLogger(__name__)


class SeuratV4Extractor(SeuratExtractor):
    """Extractor for Seurat v4/Legacy objects (uses GetAssayData API)."""

    def _extract_assay(self, assay_name: str, layer_name: str) -> tuple:
        """
        Extract sparse matrix using Seurat v4's GetAssayData function.

        Note: In v4, 'layer_name' maps to 'slot' (data, counts, scale.data)

        Returns:
            Tuple of (X_sparse, var_names, obs_names)
        """
        logger.debug(f"V4: Extracting slot '{layer_name}' from assay '{assay_name}'")
        r(f'''
            mat <- GetAssayData({self.r_var_name}, assay="{assay_name}", slot="{layer_name}")
        ''')

        X_sparse, var_names, obs_names = self._parse_dgc_matrix()
        logger.debug(f"Extracted matrix: {X_sparse.shape[0]} genes x {X_sparse.shape[1]} cells")

        return X_sparse, var_names, obs_names

    def _add_layer(self, adata: ad.AnnData, assay_name: str, layer_name: str) -> None:
        """
        Extract an additional layer and add to adata.layers.
        """
        logger.debug(f"V4: Adding slot '{layer_name}' from assay '{assay_name}'")
        try:
            r(f'''
                mat <- GetAssayData({self.r_var_name}, assay="{assay_name}", slot="{layer_name}")
            ''')
            X_sparse, _, _ = self._parse_dgc_matrix()
            adata.layers[layer_name] = X_sparse.T
            logger.info(f"Added layer '{layer_name}' to adata.layers")
        except Exception as e:
            logger.warning(f"Failed to extract slot '{layer_name}': {e}")

    def _parse_dgc_matrix(self) -> tuple:
        """
        Parse dgCMatrix components from R's 'mat' variable.

        Returns:
            Tuple of (X_sparse, var_names, obs_names)
        """
        i = np.array(r('as.integer(mat@i)'))
        p = np.array(r('as.integer(mat@p)'))
        x = np.array(r('as.numeric(mat@x)'))
        dims = np.array(r('as.integer(mat@Dim)'))

        X_sparse = sp.csc_matrix((x, i, p), shape=(dims[0], dims[1]))

        var_names = list(r('rownames(mat)'))
        obs_names = list(r('colnames(mat)'))

        return X_sparse, var_names, obs_names
