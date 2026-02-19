"""
Bacterial Genome Analyser: genome_utils.py

Utility functions for genome analysis.
Helper functions for searching, displaying, and organising gene data.
"""

from typing import Any

def print_gene_summary(gene_dict: dict[str, Any]) -> None:
    """
    Print a gene dictionary in a formatted way.
    
    Args:
        gene_dict: Dictionary containing gene information
    """
    print(f"Gene: {gene_dict["gene"]}")
    print(f"  Product: {gene_dict["product"]}")
    print(f"  Location: {gene_dict["start"]}-{gene_dict['end']}")
    print(f"  Length: {gene_dict["length"]} bp")
    print(f"  Strand: {gene_dict["strand"]}\n")

def find_gene_by_name(
    gene_data: list[dict[str, Any]], 
    gene_name: str
) -> dict[str, Any] | None:

    """
    Find a gene by its name.
    
    Args:
        gene_data: List of gene dictionaries
        gene_name: Name of gene to find
    
    Returns:
        dict: Gene dictionary if found, None otherwise
    """
    for gene in gene_data:
        if gene["gene"].lower() == gene_name.lower():
            return gene
    return None

def categorise_genes_by_size(
    gene_data: list[dict[str, Any]]
) -> dict[str, list[dict[str, Any]]]:

    """
    Categorise genes into size bins.
    
    Args:
        gene_data: List of gene dictionaries
    
    Returns:
        Dictionary mapping size category names to lists of gene dictionaries
    """

    tiny_genes = [g for g in gene_data if g["length"] < 300]
    small_genes = [g for g in gene_data if 300 <= g["length"] < 900]
    med_genes = [g for g in gene_data if 900 <= g["length"] < 2000]
    large_genes = [g for g in gene_data if 2000 <= g["length"] < 4000]
    huge_genes = [g for g in gene_data if g["length"] >= 4000]

    return {
        "tiny_genes": tiny_genes,
        "small_genes": small_genes,
        "med_genes": med_genes,
        "large_genes": large_genes,
        "huge_genes": huge_genes
    }