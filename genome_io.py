"""
Bacterial Genome Analyser: genome_io.py

Genome input/output functions for loading and reading genome files.
"""
import logging
import gzip
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

import config

logger = logging.getLogger(__name__)

def load_sequence(path: Path, seq_type: str) -> tuple[str, SeqRecord]:
    """
    Load a genome sequence from a file.

    Args:
        path: Path to the genome file
        seq_type: File format type ('genbank')

    Returns:
        Tuple containing:
            - seq: Uppercase genome sequence string
            - record: BioPython SeqRecord object containing the full genome record

    Raises:
        ValueError: If seq_type is not supported or if sequence is empty
        FileNotFoundError: If the file at path does not exist
        RuntimeError: If BioPython cannot parse the file
    """

    valid_formats = {"genbank"}
    logger.info(f"Loading {seq_type} from '{path}'")
    
    if seq_type.lower() not in valid_formats:
        logger.error(f"Unsupported file format {seq_type}")
        raise ValueError(f"Unsupported file format: {seq_type}. "
                            f"Supported formats: {valid_formats}")
    if not path.exists():
        logger.error(f"File not found:'{path}'")
        raise FileNotFoundError(f"File not found: {path}")

    try:
        if path.suffix == ".gz":
            with gzip.open(path, "rt") as handle:
                record = SeqIO.read(handle, seq_type)
        else:
            record = SeqIO.read(path, seq_type)
    
    except Exception as e:
        logger.error(f"Error parsing {seq_type} from '{path}': {e}")
        raise RuntimeError(f"Failed to load {seq_type} sequence from {path}") from e

    seq = str(record.seq).upper()
    
    if not seq:
        logger.error(f"Empty or corrupted sequence in file: '{path}'")
        raise ValueError(
            f"Genome file contains an empty sequence. "
            f"File: {path}. Check if the file is corrupted or incomplete."
        )

    logger.info(f"Successfully loaded {len(seq):,} bp sequence")
    return seq, record

def extract_gene_info(feature: SeqFeature) -> dict[str, Any]:
    """
    Extract relevant information from a CDS feature.
    
    Args:
        feature: BioPython SeqFeature object representing a CDS feature
    
    Returns:
        Dictionary containing:
            - gene: Gene name (or 'Unknown' if not annotated)
            - product: Gene product description (or 'No product' if not annotated)
            - start: Start position in base pairs
            - end: End position in base pairs
            - length: Length of the gene in base pairs
            - strand: Strand orientation ('+' or '-')
    """
    
    if not feature.location or feature.location.start is None or feature.location.end is None:
        logger.error(f"Gene feature validation failed: location={feature.location}, type={feature.type}")
        raise ValueError("Feature has missing or malformed location")

    strand = (
        "+" if feature.location.strand == 1
        else "-" if feature.location.strand == -1
        else None
    )

    return {
        "gene":feature.qualifiers.get("gene", ["Unknown"])[0],
        "product": feature.qualifiers.get("product", ["No product"])[0],
        "start": int(feature.location.start),
        "end": int(feature.location.end),
        "length": int(len(feature.location)),
        "strand": strand
    }

def print_features_summary(features: list[SeqFeature], n: int = 5) -> None:
    """
    Print a summary of the first n features in a genome record.
    
    Args:
        features: List of BioPython SeqFeature objects
        n: Number of features to print (default: 5)
    """
    
    if not features:
        logger.warning(f"No features to print - features list is empty")
        return 

    logger.debug(f"Printing summary of first {n} features")

    for idx, feature in enumerate(features[:n], start=1):
        start = int(feature.location.start)
        end = int(feature.location.end)

        if feature.type:
            print(
                f"Feature {idx}: "
                f"type={feature.type}, "
                f"location={feature.location}, "
                f"length={end - start}")

def export_gene_data(
    gene_df: pd.DataFrame, 
    genome_name: str,
    output_dir: Path = config.PROCESSED_DATA_DIR
) -> pd.DataFrame:

    """
    Export gene data and derived summaries to CSV files.
    
    Generates four CSV files:
        - ecoli_genes.csv: Full gene dataset
        - strand_summary.csv: Summary statistics grouped by strand
        - long_genes.csv: Genes longer than 3000 bp
        - ecoli_genes_categorised.csv: Full dataset with an added size category column
    
    Args:
        gene_df: DataFrame containing gene information
        output_dir: Directory to write CSV files to (default: from config)
    
    Returns:
        Copy of gene_df with an additional 'size_category' column
    """

    logger.info(f"Starting CSV export for {genome_name}")

    csv_dir = output_dir / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)
    logger.debug(f"CSV directory created: {csv_dir.relative_to(output_dir)}")

    # Main gene data
    csv_gene = csv_dir / f"{genome_name}_genes.csv"
    safe_to_csv(gene_df, csv_gene, index=False)
    logger.info(f"Exported {len(gene_df)} genes to: {csv_gene.name}")

    print(f"\nExported gene data: {csv_gene}")

    # Summary statistics by strand
    summary = (
        gene_df
        .groupby("strand")
        .agg({"length": ["count", "mean", "median", "min", "max"]})
        .round(2)
    )
    csv_summary = csv_dir / f"{genome_name}_strand_summary.csv"
    safe_to_csv(summary, csv_summary)
    logger.info(f"Exported gene summary: {csv_gene.name}")
    # Long genes only
    csv_long_genes = csv_dir / f"{genome_name}_long_genes.csv"
    safe_to_csv(
        gene_df[gene_df["length"] > 3000],
        csv_long_genes, 
        index=False
    )
    logger.info(f"Exported {len(gene_df[gene_df["length"] > 3000])} long genes")

    # Categorised genes
    gene_df_categorised = gene_df.copy()
    gene_df_categorised["size_category"] = pd.cut(
        gene_df["length"],
        bins=config.SIZE_BINS,
        labels=config.SIZE_LABELS
    )

    csv_categorised = csv_dir / f"{genome_name}_genes_categorised.csv"
    safe_to_csv(    
        gene_df_categorised,
        csv_categorised, 
        index=False
    )
    logger.info(f"Exported categorised genes to: {csv_categorised.name}")

    logger.info(f"CSV export complete: 4 files created in {csv_dir}")
    print("\nAdditional exports complete.")

    return gene_df_categorised

def safe_to_csv(gene_df: pd.DataFrame, path: Path, **kwargs) -> None:
    """
    Write a DataFrame to a CSV file with contextual I/O error handling.

    Wraps pandas.DataFrame.to_csv and re-raises any OS-level errors with
    a clear message indicating which output file failed to write.

    Args:
        df: DataFrame to write.
        path: Output CSV file path.
        **kwargs: Additional keyword arguments passed to DataFrame.to_csv.

    Raises:
        RuntimeError: If the CSV file cannot be written.
    """

    try:
        gene_df.to_csv(path, **kwargs)

    except OSError as e:
        logger.error(f"Failed to write csv file: {path}")
        raise RuntimeError(f"Failed to write csv file: {path}") from e