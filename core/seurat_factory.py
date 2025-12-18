
from seurat_v5_extractor import SeuratV5Extractor
from seurat_v4_extractor import SeuratV4Extractor


class SeuratFactory:
    """
    Creates instances of SeuratExtractor that are version compatible to the .rds file
    """

    @staticmethod
    def get_anndata(data, assay='RNA', layer='counts'):
        """
        Orchestrates selection of version-compatible extractor using version-based naming conventions
        """

        assay_info = data.get('assays', {}).get('attributes', {})
        assay_class = assay_info.get('class', [None])[0]

        if assay_class == "Assay5":
            print(f"Detected Seurat v5 object. Using v5 extractor.")
            extractor = SeuratV5Extractor()
        else:
            print(f"Detected Seurat v4/Legacy object. Using legacy extractor.")
            extractor = SeuratV4Extractor()

        return extractor.to_anndata(data, assay, layer)
