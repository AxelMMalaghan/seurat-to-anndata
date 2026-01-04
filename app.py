#!/usr/bin/env python3
"""
Command-line interface for converting Seurat objects to AnnData format.

Usage:
    python app.py input.rds -o output.h5ad
    python app.py input.rds --layers counts data --assay RNA
"""

import argparse
import logging
import sys
from pathlib import Path

from core.seurat_factory import SeuratFactory, SeuratLoadError


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert Seurat .rds files to AnnData .h5ad format",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s data.rds                          # Convert with defaults
  %(prog)s data.rds -o output.h5ad           # Specify output path
  %(prog)s data.rds --layers data counts     # Extract multiple layers
  %(prog)s data.rds --assay integrated       # Use specific assay
  %(prog)s data.rds -v                       # Verbose output
        """
    )

    parser.add_argument(
        "input",
        type=Path,
        help="Path to input .rds file containing Seurat object"
    )

    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=None,
        help="Output .h5ad file path (default: input name with .h5ad extension)"
    )

    parser.add_argument(
        "--assay",
        type=str,
        default=None,
        help="Assay to extract (default: active assay)"
    )

    parser.add_argument(
        "--layers",
        type=str,
        nargs="+",
        default=["data"],
        help="Layers to extract; first becomes X, rest go to adata.layers (default: data)"
    )

    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose (DEBUG) logging"
    )

    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress all output except errors"
    )

    return parser.parse_args()


def setup_logging(verbose: bool, quiet: bool):
    if quiet:
        level = logging.ERROR
    elif verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logging.basicConfig(
        level=level,
        format="%(levelname)s: %(message)s"
    )


def main():
    args = parse_args()
    setup_logging(args.verbose, args.quiet)

    logger = logging.getLogger(__name__)

    # Determine output path
    if args.output is None:
        output_path = args.input.with_suffix(".h5ad")
    else:
        output_path = args.output

    try:
        # Load and convert
        logger.info(f"Loading: {args.input}")
        extractor = SeuratFactory.load(str(args.input))

        logger.info(f"Converting with layers: {args.layers}")
        adata = extractor.to_anndata(
            assay_name=args.assay,
            layers=args.layers
        )

        # Save
        logger.info(f"Saving to: {output_path}")
        adata.write_h5ad(output_path)

        logger.info(f"Done! Created {adata.n_obs} cells x {adata.n_vars} genes")
        return 0

    except SeuratLoadError as e:
        logger.error(f"Failed to load Seurat object: {e}")
        return 1
    except Exception as e:
        logger.exception(f"Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
