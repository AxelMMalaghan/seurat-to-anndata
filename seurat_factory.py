import pandas as pd
import scanpy as sc
from rds2py import read_rds
from seurat_v5_extractor import SeuratV5Extractor
from seurat_v4_extractor import SeuratV4Extractor




class SeuratFactory:
    """
    Entry point for converting and Seurat RDS to AnnData
    """

    @staticmethod
    def load(file_path, assay='RNA', layer='counts'):

        r_object = read_rds(file_path)
        data = r_object['data']

        # Version detection logic
        assay_info = data.get('assays', {}).get('attributes', {})
        assay_class = assay_info.get('class', [None])[0]

        if assay_class == "Assay5":
            # replace with logging
            print(f"Detected Seurat v5 object. Using v5 extractor")
            extractor = SeuratV5Extractor()
        else:
            print(f"Detected Seurat v4/Legacy object. Using legacy extractor.")
            extractor = SeuratV4Extractor()

        return extractor.to_anndata(data, assay, layer)

