"""

Bacterial Genome Analyser: genome_visualiser.py

Genome visualisation functions for plotting.
"""

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import config
from genome_analysis import (
    split_genes_by_strand, 
    calculate_gc_window, 
    calculate_gene_density
)

def save_plot(filename: str, output_dir: Path, dpi: int = 300, display=True) -> None:
    """
    Save the current matplotlib figure to disk and display it.
    
    Args:
        filename: Name of the output file (e.g. 'plot.png')
        output_dir: Directory to save the file to
        dpi: Resolution of the saved figure in dots per inch (default: 300)
    """

    plt.tight_layout()
    filepath = output_dir / filename
    plt.savefig(filepath, dpi=dpi, bbox_inches="tight")
    print(f"Saved: {filename}")

    if display:
        plt.show()

    plt.close()

def plot_gene_length_histogram(
    gene_df: pd.DataFrame, 
    genome_name: str,
    figsize: tuple[int, int] = config.FIG_SIZE,
) -> None:
    
    """
    Plot a histogram of gene lengths with mean and median lines.
    
    Args:
        gene_df: DataFrame containing gene information with a 'length' column
        figsize: Width and height of the figure in inches (default: from config)
    """

    plt.figure(figsize=figsize)
    plt.hist(gene_df["length"], bins=50, color="steelblue", edgecolor="black", alpha=0.7)

    # Labels and title
    plt.xlabel("Gene Length (bp)", fontsize=12)  
    plt.ylabel("Number of Genes", fontsize=12)
    plt.title(f"Distribution of Gene Lengths: {genome_name}", fontweight="bold", fontsize=14)  

    # Add grid for readability
    plt.grid(axis="y", alpha=0.3)

    # Add statistics as text
    plt.axvline(gene_df["length"].mean(), color="red", linestyle="--", linewidth=2, label=f"Mean: {gene_df["length"].mean():.0f} bp" )
    plt.axvline(gene_df["length"].median(), color="green", linestyle="--", linewidth=2, label=f"Median: {gene_df["length"].median():.0f} bp")
    plt.legend()

def plot_gene_length_by_strand(
    forward_df: pd.DataFrame, 
    reverse_df: pd.DataFrame,
    genome_name: str
) -> None:
    
    """
    Plot side-by-side histograms of gene lengths for each strand.
    
    Args:
        forward_df: DataFrame of genes on the forward (+) strand
        reverse_df: DataFrame of genes on the reverse (-) strand
    """

    # Create side-by-side histograms
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(f"Gene Length Distribution by Strand: {genome_name}", fontsize=14, fontweight="bold")

    # Forward strand
    forward_lengths = forward_df["length"]
    ax1.hist(forward_lengths, bins=50, color="dodgerblue", edgecolor="black", alpha=0.7)
    ax1.axvline(forward_lengths.mean(), color="red", linestyle="--", linewidth=2)
    ax1.set_xlabel("Gene Length (bp)", fontsize=11)
    ax1.set_title(f"Forward Strand (+)\nn={len(forward_lengths)} genes", fontweight="bold", fontsize=12)
    ax1.grid(axis="y", alpha=0.3)

    # Reverse strand
    reverse_lengths = reverse_df["length"]
    ax2.hist(reverse_lengths, bins=50, color="coral", edgecolor="black", alpha=0.7)
    ax2.axvline(reverse_lengths.mean(), color="red", linestyle="--", linewidth=2)
    ax2.set_xlabel("Gene Length (bp)", fontsize=11)
    ax2.set_ylabel("Number of Genes", fontsize=11)
    ax2.set_title(f"Reverse Strand (-)\nn={len(reverse_lengths)} genes", fontweight="bold", fontsize=12)
    ax2.grid(axis="y", alpha=0.3)

def plot_gene_length_boxplot(gene_df: pd.DataFrame, genome_name) -> None:
    """
    Plot boxplots of gene lengths grouped by size category.
    
    Note:
        Modifies gene_df in place by adding a 'size_category' column.

    Adds a size_category column to gene_df based on config.SIZE_BINS,
    then plots a boxplot for each category with gene counts labelled.
    
    Args:
        gene_df: DataFrame containing gene information with a 'length' column
    """

    # Create size categories
    # Note: this modifies the original DataFrame
    gene_df["size_category"] = pd.cut(gene_df["length"],
                                bins=config.SIZE_BINS,
                                labels=["Tiny\n(<300)", "Small\n(300-900)", "Medium\n(900-2000)", "Large\n(2000-4000)", "Huge\n(>=4000)"])

    # Create box plots
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=gene_df, x="size_category", y="length", hue="size_category", palette="Set2", legend=False)

    plt.xlabel("Gene Size Category", fontsize=12)
    plt.ylabel("Gene Length (bp)", fontsize=12)
    plt.title(f"Gene Length Distribution by Size Category: {genome_name}", fontweight="bold", fontsize=14)
    plt.grid(axis="y", alpha=0.3)

    # Add count labels
    for i, category in enumerate(gene_df["size_category"].cat.categories):
        count = len(gene_df[gene_df["size_category"] == category])
        plt.text(i, plt.ylim()[1]*0.95, f"n={count}",
                ha="center", fontweight="bold", fontsize=10)

def plot_gc_content_sliding_window(
    positions: list[int], 
    gc_values: list[float],
    genome_name: str,
) -> None:
    
    """
    Plot GC content across the genome using a sliding window.
    
    Displays a line plot of GC percentage at each window position,
    with a horizontal line marking the mean GC content. The x-axis
    is formatted in megabases.
    
    Args:
        positions: List of midpoint positions for each window in base pairs
        gc_values: List of GC content percentages for each window
    """

    # Calculate GC Content
    print("\nCalculating GC content in 10kb windows...")
    print(f"Calculated {len(positions)} windows")

    # Create GC content plot
    plt.figure(figsize=(14, 6))
    plt.plot(positions, gc_values, linewidth=0.5, alpha=0.7, color="steelblue")
    plt.axhline(y=np.mean(gc_values), color="red", linestyle="--", linewidth=2, label=f"Mean GC: {np.mean(gc_values):.2f}%")

    plt.xlabel("Genomic Position (bp)",fontsize=12)
    plt.ylabel("GC Content (%)", fontsize=12) 
    plt.title(f"GC Content Variation: {genome_name} Chromosome", fontweight="bold", fontsize=14) 
    plt.grid(alpha=0.3)
    plt.legend()

    # Format x-axis to show in Mb
    ax = plt.gca()
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x/1e6:.1f}"))
    plt.xlabel("Genomic Position (Mb)", fontsize=12)

def plot_gc_distribution(
    gc_values: list[float], 
    gc_mean: float, 
    gc_std: float, 
    positions: list[int], 
    threshold_low: float, 
    threshold_high: float,
    genome_name: str
) -> None:
    
    """
    Plot the distribution of GC content across genome windows and identify unusual regions.
    
    Prints summary statistics and any regions falling outside 2 standard deviations
    from the mean, then plots a histogram with mean and threshold lines marked.
    
    Args:
        gc_values: List of GC content percentages for each window
        gc_mean: Mean GC content across all windows
        gc_std: Standard deviation of GC content across all windows
        positions: List of midpoint positions for each window in base pairs
        threshold_low: Lower GC threshold (mean - 2σ)
        threshold_high: Upper GC threshold (mean + 2σ)
    """

    print(f"\nGC Content Statistics:")
    print(f"    Mean: {gc_mean:.2f}%")
    print(f"    Std Dev: {gc_std:.2f}%")
    print(f"    Range: {min(gc_values):.2f}% - {max(gc_values):.2f}%")

    unusual_regions = []
    for i, (pos, gc) in enumerate(zip(positions, gc_values)):
        if gc > threshold_high or gc < threshold_low:
            unusual_regions.append((pos, gc))

    print(f"\nUnusual GC regions (>2σ from mean): {len(unusual_regions)}")
    if unusual_regions:
        print("\nFirst 5 unusual regions:")
        for pos, gc in unusual_regions[:5]:
            print(f"    Position: {pos:,} bp, GC: {gc:.2f}%")
    
    # GC content distribution
    plt.figure(figsize=config.FIG_SIZE)
    plt.hist(gc_values, bins=50, color="green", edgecolor="black", alpha=0.7)
    plt.axvline(gc_mean, color="red", linestyle="--", linewidth=2, label=f"Mean: {gc_mean:.2f}%")
    plt.axvline(threshold_high, color="orange", linestyle=":", linewidth=2, label=f"High threshold: {threshold_high:.2f}%")
    plt.axvline(threshold_low, color="orange", linestyle=":", linewidth=2, label=f"Low threshold: {threshold_low:.2f}%")

    plt.xlabel("GC Content (%)", fontsize=12)
    plt.ylabel("Number of Windows", fontsize=12)
    plt.title(f"Distribution of GC Content Across Genome Windows: {genome_name}", fontweight="bold", fontsize=14)
    plt.legend()
    plt.grid(axis="y", alpha=0.3)

def plot_gene_density(
    gene_positions: list[int], 
    gene_density:list[int],
    genome_name:str
) -> None:
    """
    Plot gene density across the genome.
    
    Prints summary statistics, then displays a line plot of gene counts
    per window along the chromosome. The x-axis is formatted in megabases.
    
    Args:
        gene_positions: List of midpoint positions for each window in base pairs
        gene_density: List of gene counts per window
    """

    print(f"Calculated {len(gene_positions)} windows")
    print(f"Gene density range:  {min(gene_density)} - {max(gene_density)} genes per 10kb")
    print(f"Average: {np.mean(gene_density):.2f} gene per 10kb")
    
    # Plot gene density
    plt.figure(figsize=(14, 6))
    plt.plot(gene_positions, gene_density, linewidth=0.8, alpha=0.7, color="darkgreen")
    plt.axhline(y=np.mean(gene_density), color="red", linestyle="--", linewidth=2,
                                label=f"Mean: {np.mean(gene_density):.2f} genes/10kb")

    plt.xlabel("Genomic Position (Mb)", fontsize=12)
    plt.ylabel("Number of Genes per 10kb" , fontsize=12)
    plt.title(f"Gene Density Across Chromosome: {genome_name}", fontsize=14, fontweight="bold")
    plt.grid(alpha=0.3)
    plt.legend()

    # Format x-axis
    ax = plt.gca()
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x/1e6:.1f}"))

def plot_gc_and_density_combined(
    gene_positions: list[int], 
    gene_density: list[int], 
    positions: list[int], 
    gc_values: list[float],
    genome_name: str
) -> None:
    
    """
    Plot gene density and GC content together on a dual-axis chart.
    
    Gene density is plotted on the left y-axis and GC content on the
    right y-axis, allowing direct visual comparison across the chromosome.
    The x-axis is formatted in megabases.
    
    Args:
        gene_positions: List of midpoint positions for gene density windows in base pairs
        gene_density: List of gene counts per window
        positions: List of midpoint positions for GC content windows in base pairs
        gc_values: List of GC content percentages for each window
    """

    # Create dual-axis plot
    fig, ax1 = plt.subplots(figsize=(14, 7))
    # Gene density on left axis
    color1 = "darkgreen"
    ax1.set_xlabel("Genomic Position (Mb)", fontsize=12)
    ax1.set_ylabel("Genes per 10kb", fontsize=12, color=color1)
    ax1.plot(positions, gene_density, linewidth=0.8, alpha=0.7, color=color1, label="Gene Density")
    ax1.tick_params(axis="y", labelcolor=color1)
    ax1.grid(alpha=0.3)

    # GC content on right axis
    ax2 = ax1.twinx()
    color2 = "steelblue"
    ax2.set_ylabel("GC Content (%)", color=color2, fontsize=12)
    ax2.plot(positions, gc_values, linewidth=0.5, alpha=0.5, color=color2, label="GC Content")
    ax2.tick_params(axis="y", labelcolor=color2)

    # Format x-axis
    ax1.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x/1e6:.1f}"))

    plt.title(f"Gene Density and GC Content Across Chromosome: {genome_name}", fontsize=14, fontweight="bold")

    # Combined legend
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper right")

def plot_genome_summary_figure(
    gene_df: pd.DataFrame, 
    positions: list[int], 
    gene_positions: list[int], 
    gc_values: list[float], 
    gc_mean: float, 
    gene_density: list[int],
    genome_name: str
) -> None:
    
    """
    Plot a 3-panel summary figure of the genome.
    
    Panel A shows GC content variation, Panel B shows gene density,
    and Panel C shows gene positions coloured by strand orientation.
    All panels share a common x-axis formatted in megabases.
    
    Args:
        gene_df: DataFrame containing gene information
        positions: List of midpoint positions for GC content windows in base pairs
        gene_positions: List of midpoint positions for gene density windows in base pairs
        gc_values: List of GC content percentages for each window
        gc_mean: Mean GC content across all windows
        gene_density: List of gene counts per window
    """

    # Create a 3-panel summary figure
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(16, 22), sharex=True)

    # Panel 1: GC Content
    ax1.plot(positions, gc_values, linewidth=0.5, color="steelblue", alpha=0.7)
    ax1.axhline(gc_mean, color="red", linestyle="--", linewidth=1.5)
    ax1.set_ylabel("GC Content (%)", fontsize=11)
    ax1.set_title("A. GC Content Variation", loc="left", fontsize=12, fontweight="bold")
    ax1.grid(alpha=0.3)

    # Panel 2: Gene Density
    ax2.plot(gene_positions, gene_density, linewidth=0.8, color="darkgreen", alpha=0.7)
    ax2.axhline(np.mean(gene_density), color="red", linestyle="--", linewidth=1.5)
    ax2.set_ylabel("Genes per 10kb", fontsize=11)
    ax2.set_title("B. Gene Density", loc="left", fontsize=12, fontweight="bold")
    ax2.grid(alpha=0.3)

    # Panel 3: Gene positions (scatter by strand)
    forward, reverse = split_genes_by_strand(gene_df)
    ax3.scatter(forward["start"], [1]*len(forward), alpha=0.3, s=1, c="blue", label="Forward (+)")
    ax3.scatter(reverse["start"], [-1]*len(reverse), alpha=0.3, s=1, c="red", label="Reverse (-)")
    ax3.set_yticks([-1, 1])
    ax3.set_yticklabels(["Reverse", "Forward"])
    ax3.legend(markerscale=5)
    ax3.grid(alpha=0.3)

    # Format x-axis
    ax3.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x/1e6:.1f}"))
    plt.suptitle(f"Genome Wide Analysis: {genome_name}", y=0.995, fontsize=16,  fontweight="bold")

def generate_all_visualisations(
    seq: str, 
    gene_df: pd.DataFrame,
    genome_name:str, 
    output_dir=None,
    display_plots=True,
) -> None:

    """
    Generate and save all genome visualisations.
    
    Produces 8 plots in total:
        - Gene length distribution histogram
        - Gene length by strand (side-by-side histograms)
        - Gene length boxplot by size category
        - GC content sliding window
        - GC content distribution
        - Gene density across chromosome
        - Combined GC content and gene density
        - 3-panel genome summary figure
    
    Args:
        seq: Uppercase genome sequence string
        gene_df: DataFrame containing gene information
    """

    if not display_plots:
        import matplotlib
        matplotlib.use("Agg")

    if output_dir is None:
        vis_dir = output_dir
    else:
        vis_dir = output_dir / "visualisations"
        vis_dir.mkdir(parents=True, exist_ok=True)

    # Gene length plots
    print("\nGenerating gene length visualisations...")
    forward_df, reverse_df = split_genes_by_strand(gene_df)
    plot_gene_length_histogram(gene_df, genome_name, figsize=config.FIG_SIZE)
    save_plot(f"{genome_name}_gene_length_distribution.png", vis_dir, display=display_plots)
    plot_gene_length_by_strand(forward_df, reverse_df, genome_name)
    save_plot(f"{genome_name}_gene_length_by_strand.png", vis_dir, display=display_plots)
    plot_gene_length_boxplot(gene_df, genome_name)
    save_plot(f"{genome_name}_gene_length_boxplot.png", vis_dir, display=display_plots)
    
    # GC content plots
    print("\nGenerating GC content visualisations...")
    positions, gc_values = calculate_gc_window(seq, config.WINDOW_SIZE, config.STEP_SIZE)
    gc_mean = np.mean(gc_values)
    gc_std = np.std(gc_values)
    threshold_high = gc_mean + 2 * gc_std
    threshold_low = gc_mean - 2 * gc_std
    
    plot_gc_content_sliding_window(positions, gc_values, genome_name)
    save_plot(f"{genome_name}_gc_content_sliding_window.png", vis_dir, display=display_plots)
    plot_gc_distribution(gc_values, gc_mean, gc_std, positions, threshold_low, threshold_high, genome_name)
    save_plot(f"{genome_name}_gc_content_distribution.png", vis_dir, display=display_plots)
    
    # Gene density plots
    print("\nGenerating gene density visualisations...")
    gene_positions, gene_density = calculate_gene_density(
        gene_df, len(seq), config.WINDOW_SIZE, config.STEP_SIZE)

    plot_gene_density(gene_positions, gene_density, genome_name)
    save_plot(f"{genome_name}_gene_density_across_chromosome.png", vis_dir, display=display_plots)
    plot_gc_and_density_combined(gene_positions, gene_density, positions, gc_values, genome_name)
    save_plot(f"{genome_name}_gc_and_gene_density_combined.png", vis_dir, display=display_plots)
    plot_genome_summary_figure(gene_df, positions, gene_positions, gc_values, gc_mean, gene_density, genome_name)
    save_plot(f"{genome_name}_genome_summary_figure.png", vis_dir, display=display_plots)
    
    print("\nAll visualisations complete")