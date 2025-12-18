import rds2py
import numpy as np


def summarize_structure(data, indent=0, max_depth=3):
    """
    Recursively explores the dictionary structure returned by rds2py
    to reveal slots, attributes, and data shapes.
    """
    if indent // 4 > max_depth:
        print(" " * indent + "...")
        return

    if isinstance(data, dict):
        for key, value in data.items():
            # Check for common Seurat indicators
            attr_info = ""
            if isinstance(value, dict) and 'attributes' in value:
                attr_info = f" [Attributes: {list(value['attributes'].keys())}]"

            val_type = type(value).__name__
            print(f"{' ' * indent}key: {key} ({val_type}){attr_info}")
            summarize_structure(value, indent + 4, max_depth)

    elif isinstance(data, (np.ndarray, list)):
        shape = getattr(data, 'shape', len(data))
        print(f"{' ' * indent}data: {type(data).__name__} with shape/length {shape}")


def inspect_seurat_file(path):
    print(f"--- Inspecting RDS File: {path} ---")
    # parse_rds gives us the raw dictionary representation of the R object
    obj = rds2py.parse_rds(path)

    # 1. Identify Versioning
    if 'version' in obj.get('attributes', {}):
        print(f"Detected Seurat Version Attribute: {obj['attributes']['version']}")

    # 2. Identify Assays Structure (v4 slots vs v5 layers)
    assays = obj.get('data', {}).get('assays', {})
    print(f"\nAssays Found: {list(assays.get('data', {}).keys())}")

    # 3. Recursive Summary
    print("\n--- Full Schema Summary ---")
    summarize_structure(obj)


if __name__ == "__main__":
    # Path from your app.py
    demo_path = "/home/axelm@.../KH_combined_2023-Jan-11.rds"
    try:
        inspect_seurat_file(demo_path)
    except Exception as e:
        print(f"Error reading file: {e}")