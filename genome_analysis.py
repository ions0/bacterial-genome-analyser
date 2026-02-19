"""
Bacterial Genome Analyser: genome_analysis.py

Genome analysis functions for statistical calculations.
"""

import logging
import statistics
import config
from typing import Any

import pandas as pd
from Bio.SeqIO import SeqRecord
from Bio.SeqFeature import SeqFeature

logger = logging.getLogger(__name__)

def get_gc(seq: str) -> int:
    """
    Count G and C bases in a DNA sequence.
    
    Args:
        seq: DNA sequence string (case-insensitive)
    
    Returns:
        Count of G and C bases
    """
    
    return sum(1 for s in seq.lower() if s in ("g", "c"))

def calculate_gene_stats(
    gene_data: list[dict[str, Any]], 
    seq: str
) -> dict[str, float]:

    """
    Calculate statistical metrics for genes in a genome.
    
    Args:
        gene_data: List of gene dictionaries, each containing at minimum a 'length' key
        seq: Complete genome sequence string
    
    Returns:
        Dictionary containing:
            - average_length: Mean gene length in base pairs
            - median_length: Median gene length in base pairs
            - total_coding_bp: Total number of coding base pairs
            - coding_percentage: Percentage of genome that codes for genes
    """

    logger.debug("Calculating gene statistics")

    if not gene_data:
        logger.error("Cannot calculate gene statistics: gene_data is empty")

        raise ValueError(
            "Cannot calculate gene statistics: gene_data is empty. "
            "Genome may not contain any CDS features."
    )

    # Statistics calculations
    gene_lengths = [gene["length"] for gene in gene_data]
    avg_len = statistics.mean(gene_lengths)
    median_len = statistics.median(gene_lengths)
    total_coding = sum(gene_lengths)
    coding_percent = (total_coding / len(seq)) * 100

    logger.debug(f"Gene stats calculated: avg={avg_len:.1f}bp, median={median_len:.1f}bp, "
             f"coding={coding_percent:.1f}%, total_genes={len(gene_data)}")

    return {
        "average_length": avg_len,
        "median_length": median_len,
        "total_coding_bp": total_coding,
        "coding_percentage": coding_percent
    }

def calculate_gc_window(
    seq: str, window_size: int = config.WINDOW_SIZE, 
    step: int = config.STEP_SIZE
) -> tuple[list[int], list[float]]:
    
    """
    Calculate GC content across a sequence using a sliding window.
    
    Args:
        seq: DNA sequence string (case-sensitive, expects uppercase)
        window_size: Size of the sliding window in base pairs (default: from config)
        step: Step size for sliding the window in base pairs (default: from config)
    
    Returns:
        Tuple containing:
            - positions: List of midpoint positions for each window
            - gc_values: List of GC content percentages for each window
    """

    logger.debug(f"Calculating GC windows using window size: {window_size}, and step size: {step}")

    gc_values = []
    positions = []

    if len(seq) < window_size:
        logger.error(f"Sequence length ({len(seq)} bp) is shorter than "
            f"window size ({window_size} bp). Cannot calculate GC windows.")
        
        raise ValueError("Sequence too short for GC window calculation")

    for i in range(0, len(seq) - window_size, step):
        window = seq[i:i+window_size]
        gc_count = window.count("G") + window.count("C")
        gc_percent = (gc_count / window_size) * 100

        gc_values.append(gc_percent)
        positions.append(i + window_size // 2)

    return positions, gc_values

def calculate_gene_density(
    gene_df: pd.DataFrame, 
    genome_length: int, 
    window_size: int = config.WINDOW_SIZE,
    step: int = config.STEP_SIZE
) -> tuple[list[int], list[int]]:

    """
    Calculate gene density across the genome using a sliding window.
    
    Args:
        gene_df: DataFrame containing gene information with 'start' column
        genome_length: Total length of the genome in base pairs
        window_size: Size of the sliding window in base pairs (default: from config)
        step: Step size for sliding the window in base pairs (default: from config)
    
    Returns:
        Tuple containing:
            - positions: List of midpoint positions for each window
            - density: List of gene counts per window
    """
    
    logger.debug(f"Calculating gene density with window size: {window_size}, and step size of: {step}")

    if gene_df.empty:
        logger.error("Cannot calculate gene density: gene_df is empty.")
        raise ValueError("Cannot calculate gene density: gene_df is empty.")

    density = []
    positions = []

    for i in range(0, genome_length - window_size, step):
        window_start = i
        window_end = i + window_size

        # Count genes that START in this window
        genes_in_window = len(gene_df[(gene_df["start"] >= window_start) & 
                                (gene_df["start"] < window_end)])

        density.append(genes_in_window)
        positions.append(i + window_size // 2)

    return positions, density

def split_genes_by_strand(gene_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Split genes into separate DataFrames by strand orientation.
    
    Args:
        gene_df: DataFrame containing gene information with 'strand' column
    
    Returns:
        Tuple containing:
            - forward: DataFrame of genes on the forward (+) strand
            - reverse: DataFrame of genes on the reverse (-) strand
    """

    forward = gene_df[gene_df["strand"] == "+"]
    reverse = gene_df[gene_df["strand"] == "-"]

    return forward, reverse

def find_longest_cds(record: SeqRecord) -> SeqFeature:
    """
    Find the longest coding sequence (CDS) in a genome record.
    
    Args:
        record: BioPython SeqRecord object containing genome features
    
    Returns:
        SeqFeature object representing the longest CDS feature
    
    Raises:
        ValueError: If no CDS features are found in the genome
    """

    cds_features = [f for f in record.features if f.type == "CDS"]

    if not cds_features:
        logger.error("No CDS features found in genome")
        raise ValueError("No CDS features found in genome")

    return max(cds_features, key=lambda f: len(f.location))

def calculate_strand_switches(gene_data: list[dict[str, Any]]) -> int:
    """
    Count the number of times genes switch between forward and reverse strands.
    Genes are sorted by position before counting, so input order does not matter.
    
    Args:
        gene_data: List of gene dictionaries, each containing 'start' and 'strand' keys
    
    Returns:
        Number of strand switches across the genome
    """
    
    if not gene_data:
        logger.error("Cannot calculate strand switches: gene_data is empty."
            "No genes available for strand analysis.")

        raise ValueError(
            "Cannot calculate strand switches: gene_data is empty."
            "No genes available for strand analysis."
        )

    # Sort genes by position
    sorted_genes = sorted(gene_data, key=lambda g: g["start"])

    # Count strand switches
    switches = 0
    for i in range(len(sorted_genes) - 1):
        current_strand = sorted_genes[i]["strand"]
        next_strand = sorted_genes[i+1]["strand"]
        if current_strand != next_strand:
            switches += 1
    
    return switches