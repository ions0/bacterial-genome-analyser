"""
Bacterial Genome Analyser: config.py
Configuration file for Bacterial Genome Analyser

Contains paths and analysis parameters
"""
import logging
from pathlib import Path
from datetime import datetime

# Project paths
SCRIPT_PATH = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_PATH / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
CSV_DATA_DIR = PROCESSED_DATA_DIR / "csv"
VIS_DIR = PROCESSED_DATA_DIR / "visualisations"

# Default genome file
DEFAULT_GENOME = RAW_DATA_DIR / "ecoli_k12_mg1655.gbff"

# Analysis parameters
WINDOW_SIZE = 10000
STEP_SIZE = 5000

# Gene size categories
SIZE_BINS = [0, 300, 900, 2000, 4000, 10000]
SIZE_LABELS = ["Tiny", "Small", "Medium", "Large", "Huge"]

# Visualisation settings
PLOT_DPI = 300
PLOT_STYLE = 'seaborn-v0_8-darkgrid'
FIG_SIZE = (10, 6)
GREEN = "\033[32m"
RED = "\033[31m"
RESET = "\033[0m"

PHASE_1_TIPS = [
    "Check the genome file exists at the path shown above",
    "Ensure it's a GenBank format file (.gbff, .gb or .gz)",
    "Verify the path in config.py is correct",
    "Check you have read permissions on the file"
]

PHASE_2_TIPS = [
    "Check the genome file contains CDS (coding sequence) features",
    "Verify the GenBank file isn't corrupted",
    "Make sure the file has gene annotations",
    "Some files may have no protein-coding genes"
]

PHASE_3_TIPS = [
    "Check that genes have strand information (+/-)",
    "Verify gene data was extracted correctly in Phase 2",
    "Make sure there are genes to analyse"
]

PHASE_4_TIPS = [
    "Check you have write permissions in the output directory",
    "Verify there's enough disk space",
    "Make sure the data directory exists",
    "Check that gene data isn't empty"
]

PHASE_5_TIPS = [
    "This is optional - your data is already saved",
    "Check matplotlib is installed correctly",
    "Verify you have write permissions for visualisations folder",
    "If running on a server, you may need a display backend"
]

def setup_output_directories() -> None:
    """Create all required output directories"""

    directories = [
        VIS_DIR,
        CSV_DATA_DIR,
        PROCESSED_DATA_DIR,
    ]
    
    for directory in directories:
        try:
            directory.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            raise RuntimeError(f"Failed to create directory: {directory}") from e
    
    print("Output directories initialised")

def setup_logging(log_dir: Path) -> None:
    
    log_dir.mkdir(parents=True, exist_ok=True)    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"bacterial_genome_analyser_{timestamp}.log"

    LOG_FORMAT = "[%(asctime)s] [%(name)s] [%(levelname)s] [%(filename)s:%(lineno)d] [%(message)s]"
    formatter = logging.Formatter(LOG_FORMAT)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    # Silence noisy third-party libraries
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('PIL').setLevel(logging.WARNING)
    logging.getLogger('urllib3').setLevel(logging.WARNING)