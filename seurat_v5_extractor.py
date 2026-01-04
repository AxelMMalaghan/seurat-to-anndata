import numpy as np
import scipy.sparse as sp
from rpy2.robjects import r

from core.seurat_extractor import SeuratExtractor


class SeuratV5Extractor(SeuratExtractor):
    """Extractor for Seurat v5 objects (uses LayerData API)."""

    def _extract_assay(self, assay_name: str, layer_name: str) -> tuple:
        """
        Extract sparse matrix using Seurat v5's LayerData function.

        Returns:
            Tuple of (X_sparse, var_names, obs_names)
        """
        r(f'''
            mat <- LayerData({self.r_var_name}, assay="{assay_name}", layer="{layer_name}")
        ''')

        # Extract dgCMatrix components directly to avoid rpy2 conversion issues
        i = np.array(r('as.integer(mat@i)'))
        p = np.array(r('as.integer(mat@p)'))
        x = np.array(r('as.numeric(mat@x)'))
        dims = np.array(r('as.integer(mat@Dim)'))

        X_sparse = sp.csc_matrix((x, i, p), shape=(dims[0], dims[1]))

        var_names = list(r('rownames(mat)'))
        obs_names = list(r('colnames(mat)'))

        return X_sparse, var_names, obs_names
