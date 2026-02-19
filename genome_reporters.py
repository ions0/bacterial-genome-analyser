"""
    Bacterial Genome Analyser: genome_reporters.py

    Console output functions for displaying genome analysis results.
    Formats and prints summaries for each analysis phase including basic
    genome info, gene statistics, strand distribution, and DataFrame summaries.
"""
import statistics
from typing import Any

import pprint
import pandas as pd
from collections import Counter
from Bio.SeqIO import SeqRecord
from Bio.SeqFeature import SeqFeature

from genome_io import print_features_summary
from genome_utils import print_gene_summary

def print_basic_genome_info(
    record: SeqRecord, 
    seq: str, 
    gc_percent: float, 
    longest_cds: SeqFeature
) -> None:

    """
    Print a high-level summary of the genome.
    
    Displays genome ID, description, length, feature types, GC content,
    and details of the longest CDS.
    
    Args:
        record: BioPython SeqRecord object containing the genome record
        seq: Uppercase genome sequence string
        gc_percent: GC content of the genome as a percentage
        longest_cds: SeqFeature object representing the longest CDS in the genome
    """
    
    print("ID: " + record.id)
    print(f"DESC: {record.description}")
    print(f"Length: {len(seq)}")
    print(f"Features: {record.features[:10]}")

    print_features_summary(record.features, n=5)

    feature_types = [feature.type for feature in record.features]
    feature_counts = Counter(feature_types)
    pprint.pprint(feature_counts)

    print(f"GC: {gc_percent:.2f}%")

    print(f"\nLongest Gene:")
    print(f"    Length: {len(longest_cds.location)} bp")
    print(f"    Location: {longest_cds.location}")
    # Optional: print gene name if available
    if "gene" in longest_cds.qualifiers:
        print(f"    Gene: {longest_cds.qualifiers["gene"][0]}")

def print_gene_analysis_results(
    gene_data: list[dict[str, Any]], 
    genome_stats: dict[str, Any], 
    gene_sizes: dict[str, Any], 
    cds_features: list[SeqFeature]
) -> None:

    """
    Print a detailed summary of gene analysis results.
    
    Displays the first 5 genes, longest and shortest genes, genome-wide
    statistics, size distribution.
    
    Args:
        gene_data: List of gene dictionaries containing gene information
        genome_stats: Dictionary of genome statistics as returned by calculate_gene_stats()
        gene_sizes: Dictionary of genes categorised by size as returned by categorise_genes_by_size()
        cds_features: List of BioPython SeqFeature objects representing CDS features
    """

    print(f"\nTotal CDS features: {len(cds_features)}")
    pprint.pprint(cds_features[0])

    gene_name = cds_features[0].qualifiers.get("gene", ["Unknown"])[0]
    print(gene_name)

    for gene in gene_data[:5]:
        print_gene_summary(gene)

    longest_gene = max(gene_data, key=lambda g: g["length"])
    shortest_gene = min(gene_data, key=lambda g: g["length"])

    print(f"\n===== Longest Gene =====\n")
    print_gene_summary(longest_gene)
    print(f"\n===== Shortest Gene =====\n")
    print_gene_summary(shortest_gene)
    print(f"\n===== Average Gene Length =====\n")
    print(genome_stats["average_length"])
    print(f"\n===== Median Gene Length =====\n")
    print(genome_stats["median_length"])
    print(f"\n===== Total Coding Bases =====\n")
    print(genome_stats["total_coding_bp"])
    print(f"\n===== Percent of Genome that codes for genes =====\n")
    print(genome_stats["coding_percentage"])
    print(f"\n===== Gene Size Distribution =====\n")
    print(f"Tiny Genes (<300 bp): {len(gene_sizes["tiny_genes"])}")
    print(f"Small genes (300-900 bp): {len(gene_sizes["small_genes"])}")
    print(f"Medium Genes (900-2000 bp): {len(gene_sizes["med_genes"])}")
    print(f"Large Genes (2000-4000 bp): {len(gene_sizes["large_genes"])}")
    print(f"Huge Genes ( >= 4000 bp): {len(gene_sizes["huge_genes"])}")

def print_strand_analysis_results(
    forward_genes: list[dict[str, Any]], 
    reverse_genes: list[dict[str, Any]], 
    gene_data: list[dict[str, Any]], 
    switches: int
) -> None:

    """
    Print a detailed summary of strand distribution and analysis.
    
    Displays gene counts, percentages, and average lengths by strand,
    the longest gene on each strand, and strand switching statistics.
    
    Args:
        forward_genes: List of gene dictionaries on the forward (+) strand
        reverse_genes: List of gene dictionaries on the reverse (-) strand
        gene_data: List of all gene dictionaries
        switches: Total number of strand switches across the genome
    """

    print(f"\n===== Strand Distribution =====")
    print(f"    Forward strand (+): {len(forward_genes)} genes")
    print(f"    Reverse strand (-): {len(reverse_genes)} genes")

    if len(reverse_genes) > 0:
        ratio = len(forward_genes) / len(reverse_genes)
        print(f"    Ratio: {ratio:.2f}")
    else:
        print(f"    Ratio: All genes on forward strand.")

    # Calculate percentages
    forward_pc = (len(forward_genes) / len(gene_data)) * 100
    reverse_pc = (len(reverse_genes) / len(gene_data)) * 100

    print(f"\nPercentages:")
    print(f"    Forward: {forward_pc:.2f}%")
    print(f"    Reverse: {reverse_pc:.2f}%")

    # Average length by strand
    if forward_genes:
        forward_avg = statistics.mean([g["length"] for g in forward_genes])
        print(f"    Forward strand: {forward_avg:.1f} bp")
    else:
        print(f"    Forward strand: No genes")

    if reverse_genes:
        reverse_avg = statistics.mean([g["length"] for g in reverse_genes])
        print(f"    Reverse strand: {reverse_avg:.1f} bp")
    else:
        print(f"    Reverse strand: No genes")    

    # Find longest on each strand
    if forward_genes:
        longest_forward = max(forward_genes, key=lambda g: g["length"])
        print(f"\nLongest gene on forward strand:")
        print(f"    {longest_forward["gene"]}: {longest_forward["length"]} bp")
    else:
        print(f"\nLongest gene on forward strand: No genes on forward strand")

    if reverse_genes:
        longest_reverse = max(reverse_genes, key=lambda g: g["length"])
        print(f"\nLongest gene on reverse strand:")
        print(f"    {longest_reverse["gene"]}: {longest_reverse["length"]} bp")
    else:
        print(f"\nLongest gene on reverse strand: No genes on reverse strand")
         
    print(f"\nStrand Switching:")
    print(f"    Total strand switches: {switches}")

    if switches > 0:
        avg_before_switch = len(gene_data) / switches
        print(f"    Average genes before switch: {avg_before_switch:.1f}")
    else:
        print(f"    Average genes before switch: No strand switches detected")

def print_dataframe_summary(gene_df: pd.DataFrame) -> None:
    """
    Print a comprehensive summary of the gene DataFrame.
    
    Displays shape, first 5 rows, column info, length statistics,
    top 10 longest and shortest genes, strand distribution, genes
    longer than 3000 bp, ribosomal genes, and genes in the first 100kb.
    
    Args:
        gene_df: DataFrame containing gene information
    """
    
    print("\nDataFrame created")
    print(f"Shape: {gene_df.shape[0]} rows x {gene_df.shape[1]} columns")
    print("\nFirst 5 rows:")
    print(gene_df.head())

    print("\nColumn names:")
    print(gene_df.columns.tolist())

    print("\nDataFrame info:")
    print(gene_df.info())

    # Basic statistics
    print("\nLength statistics:")
    print(gene_df["length"].describe())

    # Top 10 longest genes
    print("\nTop 10 Longest Genes:")
    print(gene_df.nlargest(10, "length")[["gene", "product", "length", "strand"]])

    # Top 10 shortest genes
    print("\nTop 10 Shortest Genes:")
    print(gene_df.nsmallest(10, "length")[["gene", "product", "length", "strand"]])

    # Genes by strand
    print("\nGenes by strand:")
    print(gene_df["strand"].value_counts())

    # Average length by strand
    print("\nAverage length by strand:")
    print(gene_df.groupby("strand")["length"].mean())

    # Find long genes (> 3000 bp)
    long_genes = gene_df[gene_df["length"] > 3000]
    print(f"\nGenes longer than 3000 bp: {len(long_genes)}")
    print(long_genes[["gene", "product", "length"]])

    # Find genes with specific keywords in product
    print("\nGenes related to 'ribosom':")
    ribosomal = gene_df[gene_df["product"].str.contains("ribosom", case=False, na=False)]
    print(f"Found {len(ribosomal)} ribosomal genes")
    print(ribosomal[["gene", "product", "length"]].head(10))

    # Find genes in a specific region (e.g., first 100kb)
    region_genes = gene_df[gene_df["start"] < 100000]
    print(f"\nGenes in first 100kb: {len(region_genes)}")