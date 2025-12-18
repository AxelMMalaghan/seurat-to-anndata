import rds2py
import numpy as np
import os


def summarize_structure(data, indent=0, max_depth=3):
    if indent // 4 > max_depth:
        return

    if isinstance(data, dict):
        for key, value in data.items():
            attr_info = ""
            if isinstance(value, dict) and 'attributes' in value:
                # This identifies the Seurat Assay version (v4 vs v5)
                class_attr = value['attributes'].get('class', [None])[0]
                attr_info = f" [Class: {class_attr}]"

            print(f"{' ' * indent}key: {key} ({type(value).__name__}){attr_info}")
            summarize_structure(value, indent + 4, max_depth)

    elif isinstance(data, (np.ndarray, list)):
        shape = getattr(data, 'shape', len(data))
        print(f"{' ' * indent}data: {type(data).__name__} with shape/length {shape}")


def inspect_demo_file(path):
    if not os.path.exists(path):
        print(f"Error: The file does not exist at: {path}")
        return

    print(f"--- Inspecting: {path} ---")
    try:
        # Using parse_rds to read the file directly into a Python dictionary
        obj = rds2py.parse_rds(path)

        # Check for Seurat version in the object attributes
        version = "Unknown"
        if 'attributes' in obj and 'version' in obj['attributes']:
            version = obj['attributes']['version']
        print(f"Reported Seurat Version: {version}")

        # Scan the assays to find the data structure
        print("\n--- Assay Structure Scan ---")
        if 'data' in obj and 'assays' in obj['data']:
            summarize_structure(obj['data']['assays'], max_depth=4)
        else:
            print("No assays found in the expected location.")

    except Exception as e:
        print(f"RDS Parser Error: {e}")


if __name__ == "__main__":
    # The absolute path we confirmed earlier
    target_path = "/home/axelm@malaghan.org.nz/seurat-to-anndata/data/rds/KH_combined_2023-Jan-11.rds"
    inspect_demo_file(target_path)

