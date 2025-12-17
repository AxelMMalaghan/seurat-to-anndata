import pandas as pd
import scanpy as sc
from rds2py import read_rds

from dataset_loader import SeuratLoader
from seurat_v5_extractor import SeuratV5Extractor
from seurat_v4_extractor import SeuratV4Extractor




class SeuratFactory:
    """
    Orchestrates selection of version-compatible extractor
    """

    @staticmethod
    def convert(file_path, assay='RNA', layer='counts'):
        # Use the new Loader class
        loader = SeuratLoader(file_path)
        data = loader.load_data()

        # Version detection logic
        assay_info = data.get('assays', {}).get('attributes', {})
        assay_class = assay_info.get('class', [None])[0]

        if assay_class == "Assay5":
            print(f"Detected Seurat v5 object. Using v5 extractor.")
            extractor = SeuratV5Extractor()
        else:
            print(f"Detected Seurat v4/Legacy object. Using legacy extractor.")
            extractor = SeuratV4Extractor()

        return extractor.to_anndata(data, assay, layer)

