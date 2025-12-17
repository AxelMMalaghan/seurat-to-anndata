from rds2py import read_rds

class SeuratLoader:
    """
    Handles the loading of a Seurat RDS file
    """
    def __init__(self, seurat_file_path):
        self.seurat_file_path = seurat_file_path

    def load_data(self):

        try:
            r_object = read_rds(self.seurat_file_path)
            # Ensure the structure matches expectations
            if 'data' not in r_object:
                raise ValueError("Invalid Seurat RDS file")
            return r_object['data']
        except Exception as e:
            raise IOError("Seurat RDS file could not be loaded")

