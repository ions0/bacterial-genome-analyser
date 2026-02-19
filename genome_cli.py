"""
Bacterial Genome Analyser: genome_cli.py

Command-line interface module for argument parsing and validation.
"""
import logging
import sys
import argparse
from datetime import datetime
from pathlib import Path

logger = logging.getLogger(__name__)

def parse_arguments():
    """
    Parse and validate command-line arguments.
    
    Returns:
        argparse.Namespace: Parsed arguments containing:
            - genome: Path to genome file
            - output: Output directory path
    """
    
    parser = argparse.ArgumentParser(
        description="Bacterial Genome Analyser - Comprehensive analysis of bacterial genomes",
        epilog="Example: python bacterial_genome_analyser.py --genome data/raw/ecoli.gbff"
    )

    parser.add_argument(
        "--genome",
        type=str,
        required=True,
        help="Path to genome file (.gbff, .gb)"
    )

    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output directory for results (default: data/processed)"
    )

    parser.add_argument(
        "--version",
        action="version",
        version="Bacterial Genome Analyser v1.0.1"
    )

    parser.add_argument(
        "--no-display",
        action="store_true",
        help="Save plots without displaying them"
    )

    return parser.parse_args()

def validate_arguments(args):
    """
    Validate parsed arguments before running analysis.
    
    Args:
        args: Parsed arguments from parse_arguments()
    
    Raises:
        ValueError: If genome file doesn't exist or has wrong extension
        
    Returns:
        dict: Validated configuration with Path objects
    """

    genome_path = Path(args.genome)
    if not genome_path.exists():
        raise ValueError(f"Genome file not found: {genome_path}")

    valid_extensions = {".gbff", ".gb"}

    if genome_path.suffix.lower() not in valid_extensions:
        raise ValueError(f"Invalid file extension. Expected: {valid_extensions}")

    file_format = "genbank"

    # Output argument checking
    if args.output is None:
        output_dir = None
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        genome_name = genome_path.stem

        base_output = Path(args.output)
        output_dir = base_output / f"{genome_name}_{timestamp}"

    # Return clean config
    return {
        "genome_path": genome_path,
        "file_format": file_format,
        "output_dir": output_dir,
        "display_plots": not args.no_display
    }


