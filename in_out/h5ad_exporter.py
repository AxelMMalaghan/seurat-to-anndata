import os
from pathlib import Path

class H5ADExporter:
    """
    Handles the saving of AnnData objects to the local file system - for now - in .h5ad format
    """
    def __init__(self, base_dir="data/h5ad"):
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(parents=True, exist_ok=True)

    def export(self, adata, filename):
        """
        Saves the provided AnnData object to the configured directory
        """

        if not filename.endswith(".h5ad"):
            filename += ".h5ad"

        save_path = self.base_dir / filename

        print(f"Exporting AnnData object to {save_path}")
        adata.write_h5ad(save_path)
        return str(save_path)
