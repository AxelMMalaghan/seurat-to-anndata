from rpy2.robjects import r
from rpy2.robjects.packages import importr

from seurat_v5_extractor import SeuratV5Extractor
from seurat_v4_extractor import SeuratV4Extractor


class SeuratFactory:
    """
    Factory that loads a Seurat RDS file and returns an appropriate extractor
    based on the detected Seurat version.
    """

    @staticmethod
    def load(rds_path: str, r_var_name: str = "seurat_obj"):
        """
        Load an RDS file and return the appropriate version-specific extractor.

        Args:
            rds_path: Path to the .rds file containing the Seurat object
            r_var_name: Variable name to use in R's global environment

        Returns:
            SeuratExtractor subclass instance (V4 or V5)
        """
        # Ensure Seurat is loaded
        importr('Seurat')

        # Load the object into R
        r(f'{r_var_name} <- readRDS("{rds_path}")')

        # Detect version by checking assay class
        active_assay = r(f'DefaultAssay({r_var_name})')[0]
        assay_class = r(f'class({r_var_name}@assays${active_assay})')[0]

        if assay_class == "Assay5":
            print(f"Detected Seurat v5 object (Assay5). Using V5 extractor.")
            return SeuratV5Extractor(r_var_name)
        else:
            print(f"Detected Seurat v4/Legacy object ({assay_class}). Using V4 extractor.")
            return SeuratV4Extractor(r_var_name)
