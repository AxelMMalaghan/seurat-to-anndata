import logging
import os

from rpy2.robjects import r
from rpy2.robjects.packages import importr

from seurat_v5_extractor import SeuratV5Extractor
from seurat_v4_extractor import SeuratV4Extractor

logger = logging.getLogger(__name__)


class SeuratLoadError(Exception):
    """Raised when a Seurat object cannot be loaded or validated."""
    pass


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

        Raises:
            SeuratLoadError: If file doesn't exist, isn't readable, or isn't a valid Seurat object
        """
        # Validate file exists
        if not os.path.exists(rds_path):
            raise SeuratLoadError(f"File not found: {rds_path}")

        if not os.path.isfile(rds_path):
            raise SeuratLoadError(f"Path is not a file: {rds_path}")

        if not rds_path.lower().endswith('.rds'):
            logger.warning(f"File does not have .rds extension: {rds_path}")

        logger.info(f"Loading Seurat object from: {rds_path}")

        # Ensure Seurat is loaded
        try:
            importr('Seurat')
        except Exception as e:
            raise SeuratLoadError(f"Failed to load Seurat R package: {e}")

        # Load the object into R
        try:
            r(f'{r_var_name} <- readRDS("{rds_path}")')
        except Exception as e:
            raise SeuratLoadError(f"Failed to read RDS file: {e}")

        # Validate it's a Seurat object
        obj_class = r(f'class({r_var_name})')[0]
        if obj_class != "Seurat":
            raise SeuratLoadError(
                f"Object is not a Seurat object (got class: {obj_class})"
            )

        # Detect version by checking assay class
        active_assay = r(f'DefaultAssay({r_var_name})')[0]
        assay_class = r(f'class({r_var_name}@assays${active_assay})')[0]

        if assay_class == "Assay5":
            logger.info(f"Detected Seurat v5 object (Assay5). Using V5 extractor.")
            return SeuratV5Extractor(r_var_name)
        else:
            logger.info(f"Detected Seurat v4/Legacy object ({assay_class}). Using V4 extractor.")
            return SeuratV4Extractor(r_var_name)
