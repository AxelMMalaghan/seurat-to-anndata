from abc import ABC, abstractmethod

import pandas as pd
import scanpy as sc


class SeuratExtractor(ABC):
    """
    Base interface for all Seurat extractors (based on version)
    """

    def to_anndata(self, data, assay_name, layer_name):
        """
        Orchestrates the decomposition of the Seurat object by delegating to:

        * _extract_obs
        * _extract_assay

        and then recomposes the AnnData object and adds embeddings to the AnnData object (if they exist) using:

        * _add_embeddings
        """

        obs = self._extract_obs(data)
        var_names, X = self._extract_assay(data, assay_name, layer_name)

        # Build object
        adata = sc.AnnData(X=X, obs=obs)
        adata.var_names = var_names

        # Add embeddings
        adata = self._add_embeddings(adata, data)
        return adata


    @abstractmethod
    def _extract_assay(self, data, assay_name, layer_name):
        """
        Must be implemented in subclasses - this is where Seurat versions differ (from what I can tell)
        """

    def _extract_obs(self, data):

        meta = data['meta.data']
        df = pd.DataFrame(meta['data'])
        df.index = meta['attributes']['row.names']

        for col in df.columns:
            col_raw = meta['data'][col]
            if isinstance(col_raw, dict) and 'attributes' in col_raw:
                levels = col_raw['attributes'].get['levels']
                if levels:
                    df[col] = [levels[i-1] for i in col_raw['data']]

        return df


    def _add_embeddings(self, adata, data):
        """
        Adds the embeddings to the AnnData object - .obsm
        """
        if 'reductions' in data:
            for name, red, in data['reductions']['data'].items():
                adata.obsm[f"X_{name.lower}"] = red['data']['cell.embeddings']

        return adata
