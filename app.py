import scanpy as sc
import rds2py

from core.seurat_factory import SeuratFactory
from in_out.h5ad_exporter import H5ADExporter
from in_out.seurat_loader import SeuratLoader


def main(input_path, output_name):
    # 1. Load the raw R object
    loader = SeuratLoader(input_path)
    data = loader.load_data()

    # 2. Convert to AnnData
    adata = SeuratFactory.get_anndata(data)

    # 3. Export to disk
    exporter = H5ADExporter()
    saved_path = exporter.export(adata, output_name)

    print(f"Successfully converted and saved to: {saved_path}")


if __name__ == "__main__":
    main(input_path="/home/axelm@malaghan.org.nz/seurat-to-anndata/data/rds/KH_combined_2023-Jan-11.rds", output_name="processed_cells")