"""
Bacterial Genome Analyser: bacterial_genome_analyser.py
Author: Jared Cambridge
Date: January 10, 2026
Updated: February 19, 2026

Version: 1.0.1
"""

import sys
import logging
from pathlib import Path
from datetime import datetime

import pandas as pd
from Bio import SeqIO

import config
import genome_cli
from config import (
    setup_output_directories,
    setup_logging,
    PHASE_1_TIPS,
    PHASE_2_TIPS,
    PHASE_3_TIPS,
    PHASE_4_TIPS,
    PHASE_5_TIPS
)
from genome_io import load_sequence, extract_gene_info, export_gene_data
from genome_analysis import (
    get_gc, 
    calculate_gene_stats, 
    find_longest_cds, 
    calculate_strand_switches 
)  
from genome_utils import categorise_genes_by_size
from genome_visualiser import generate_all_visualisations
from genome_reporters import (
    print_basic_genome_info, 
    print_gene_analysis_results, 
    print_strand_analysis_results, 
    print_dataframe_summary
)

logger = logging.getLogger(__name__)

def show_error(phase_name, error, suggestions):
    """Show an error message"""
    print(f"\n{config.RED}{"="*50}\nERROR IN: {phase_name}\n{"="*50}{config.RESET}")
    print(f"What happened: {error}")
    print("\nWhat to try:")
    for i, suggestion in enumerate(suggestions, 1):
        print(f"  {i}. {suggestion}")
    print("="*50)

def main(genome_path=None, file_format="genbank", output_dir=None, display_plots=True):

    logger.info(f"Starting bacterial genome analysis")
    logger.info(f"Genome: {genome_path.name}, Format: {file_format}")
    logger.debug(f"Output directory: {output_dir}, Display plots: {display_plots}")

    # Set defaults if not provided
    if genome_path is None:
        logger.info(f"Using default genome: {config.DEFAULT_GENOME}")
        genome_path = config.DEFAULT_GENOME
    
    genome_name = genome_path.stem

    if output_dir is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_dir = config.PROCESSED_DATA_DIR / f"{genome_name}_{timestamp}"

    genome_title = genome_name.replace(" ", "_").title()
    
    # ==================== PHASE 1: LOAD AND ANALYSE BASIC INFO ==========================
    logger.info("Phase 1: Load & Analyse Genome")
    
    print(
        f"\n\n{config.GREEN}{"="*50}\nPhase 1: Load & Analyse Genome\n{"="*50}{config.RESET}")

    try:
        seq, record = load_sequence(genome_path, file_format)
        gc_count = get_gc(seq)
        gc_percent = gc_count / len(seq) * 100
        longest_cds = find_longest_cds(record)
        print_basic_genome_info(record, seq, gc_percent, longest_cds)

        logger.info("Phase 1 Complete")

    except Exception as e:
        logger.error(f"Phase 1 failed: {e}", exc_info=True)
        show_error("Phase 1", e, PHASE_1_TIPS)
        return

    # ==================== PHASE 2: GENE ANALYSIS ========================================
    logger.info("Phase 2: Gene Analysis")

    print(
        f"{config.GREEN}{"="*50}\nPhase 2: Gene Analysis\n{"="*50}{config.RESET}")

    try:
        cds_features = [f for f in record.features if f.type == "CDS"]
        logger.info(f"Found {len(cds_features)} CDS features")
        gene_data = [extract_gene_info(f) for f in cds_features]
        genome_stats = calculate_gene_stats(gene_data, seq)
        gene_sizes = categorise_genes_by_size(gene_data)
        print_gene_analysis_results(gene_data, genome_stats, gene_sizes, cds_features)

        logger.info("Phase 2 Complete")

    except Exception as e:
        logger.error(f"Phase 2 failed: {e}", exc_info=True)
        show_error("Phase 2", e, PHASE_2_TIPS)
        return

    # ==================== PHASE 3: STRAND ANALYSIS ======================================
    logger.info("Phase 3: Strand Analysis")

    print(
        f"{config.GREEN}{"="*50}\nPhase 3: Strand Analysis\n{"="*50}{config.RESET}")

    try:
        # Count genes on each strand
        forward_genes = [g for g in gene_data if g["strand"] == "+"]
        reverse_genes = [g for g in gene_data if g["strand"] == "-"]
        # Count switches
        switches = calculate_strand_switches(gene_data)
        logger.info(f"Strand distribution: {len(forward_genes)} forward, {len(reverse_genes)} reverse")
        print_strand_analysis_results(forward_genes, reverse_genes, gene_data, switches)
    
        logger.info("Phase 3 Complete")

    except Exception as e:
        logger.error(f"Phase 3 failed: {e}", exc_info=True)
        show_error("Phase 3", e, PHASE_3_TIPS)
        return
    
    # ================ PHASE 4: DATA ORGANISATION ========================================
    logger.info("Phase 4: Data Organisation")

    print(
        f"{config.GREEN}{"="*50}\nPhase 4: Data Organisation\n{"="*50}{config.RESET}")

    try:
        # Create DataFrame from gene_data
        gene_df = pd.DataFrame(gene_data)
        print_dataframe_summary(gene_df) 
        export_gene_data(gene_df, genome_title, output_dir=output_dir)

        logger.info("Phase 4 Complete")

    except Exception as e:
        logger.error(f"Phase 4 failed: {e}", exc_info=True)
        show_error("Phase 4", e, PHASE_4_TIPS)
        return

    # ==================== PHASE 5: VISUALISATION ========================================
    logger.info("Phase 5: Visualisation")

    print(
        f"{config.GREEN}{"="*50}\nPhase 5: Generating Visualisations\n{"="*50}{config.RESET}")
    
    try: 
        generate_all_visualisations(seq, gene_df, genome_title, output_dir, display_plots)
        
        logger.info("Phase 5 Complete")

    except Exception as e:
        logger.error(f"Phase 5 failed: {e}", exc_info=True)
        show_error("Phase 5", e, PHASE_5_TIPS)

    logger.info("Genome Analysis Complete")
    print(
        f"\n\n{config.GREEN}{"="*50}\nAnalysis Complete\n{"="*50}{config.RESET}\n\n")

if __name__ == "__main__":
    try:
        args = genome_cli.parse_arguments()
        config_dict = genome_cli.validate_arguments(args)

        # Extract genome info
        genome_name = Path(args.genome).stem.replace(" ", "_")
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        if config_dict["output_dir"] is None:
            output_dir = config.PROCESSED_DATA_DIR / f"{genome_name}_{timestamp}"  
        else:
            output_dir = config_dict["output_dir"]    

        setup_logging(output_dir / "logs")

        logger.info("Bacterial Genome Analyser v1.0.1 started")
        logger.debug(f"Command line args: {sys.argv}")
        logger.info(f"Output Directory: {output_dir}")

        main(
            genome_path=config_dict["genome_path"],
            file_format="genbank",
            output_dir=output_dir,
            display_plots=config_dict["display_plots"]
        )

    except ValueError as e:
        print(f"\n{config.RED}{"="*50}\nERROR: Invalid Arguments\n{"="*50}{config.RESET}")
        print(f"Error: {e}")
        print("\nRun with --help for usage information")
        print("="*50)
        sys.exit(1)