from in_out.seurat_loader import SeuratLoader

def main(input_path, output_name):
    loader = SeuratLoader()
    try:
        adata = loader.load_data(input_path)
        adata.write_h5ad(f"{output_name}.h5ad")
        print(f"Saved to {output_name}.h5ad")
    except Exception as e:
        print(f"Critical Error: {e}")

if __name__ == "__main__":
    main(input_path="/home/axelm@.../KH_combined_2023-Jan-11.rds",
         output_name="processed_cells")