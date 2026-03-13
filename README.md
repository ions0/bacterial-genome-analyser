# Bacterial Genome Analyser

A command-line tool for analysing bacterial genomes in GenBank format. Performs GC content analysis, gene statistics, strand distribution analysis, and generates visualisations - all exported to CSV and image files.

Built as a self-directed learning project prior to starting a bioinformatics degree.

---

## Features

- Loads GenBank (`.gbff`, `.gb`) genome files
- Calculates GC content across the full genome and in sliding windows
- Extracts and analyses all CDS (coding sequence) features
- Reports gene size distribution, strand bias, and strand switching frequency
- Calculates gene density across the genome
- Exports gene data and summary statistics to CSV
- Generates visualisations saved as high-resolution PNG files
- Structured logging to file for each run

---

## Project Structure

```
bacterial_genome_analyser/
├── bacterial_genome_analyser.py   # Main entry point
├── config.py                      # Paths and analysis parameters
├── genome_cli.py                  # Argument parsing and validation
├── genome_io.py                   # File loading and CSV export
├── genome_analysis.py             # Statistical analysis functions
├── genome_utils.py                # Helper and utility functions
├── genome_reporters.py            # Console output formatting
├── genome_visualiser.py           # Plot generation
├── data/
│   ├── raw/                       # Place genome files here
│   └── processed/                 # Output files written here
├── requirements.txt
└── README.md
```

---

## Requirements

- Python 3.11+
- See `requirements.txt` for dependencies

---

## Installation

```bash
git clone https://github.com/ions0/bacterial_genome_analyser.git
cd bacterial_genome_analyser
pip install -r requirements.txt
```

---

## Usage

```bash
python bacterial_genome_analyser.py --genome data/raw/ecoli_k12_mg1655.gbff
```

### Options

| Flag | Description |
|------|-------------|
| `--genome` | Path to genome file (required) |
| `--output` | Output directory (default: `data/processed/`) |
| `--no-display` | Save plots without displaying them |
| `--version` | Show version number |
| `--help` | Show help message |

### Example

```bash
# Basic run with auto-detected format
python bacterial_genome_analyser.py --genome data/raw/ecoli_k12_mg1655.gbff

# Specify output directory, suppress plot display
python bacterial_genome_analyser.py --genome data/raw/ecoli_k12_mg1655.gbff --output results/ --no-display
```

---

### Analysis Phases

The tool runs through five distinct phases:

1. **Phase 1: Load & Analyse Genome**
   - Loads sequence data
   - Calculates basic metrics (length, GC content)
   - Identifies longest CDS

2. **Phase 2: Gene Analysis**
   - Extracts all CDS features
   - Computes gene statistics
   - Categorizes genes by size

3. **Phase 3: Strand Analysis**
   - Analyzes strand distribution
   - Calculates strand switching patterns
   - Identifies longest genes per strand

4. **Phase 4: Data Organisation**
   - Creates structured DataFrames
   - Exports data to CSV files

5. **Phase 5: Visualisation**
   - Generates plots and figures
   - Saves to visualisations directory

---

## Output

Each run creates a timestamped folder inside the output directory containing:

```
ecoli_k12_mg1655_20260210_143022/
├── csv/
│   ├── Ecoli_K12_Mg1655_genes.csv
│   ├── Ecoli_K12_Mg1655_strand_summary.csv
│   ├── Ecoli_K12_Mg1655_long_genes.csv
│   └── Ecoli_K12_Mg1655_genes_categorised.csv
├── visualisations/
│   └── *.png
└── logs/
    └── bacterial_genome_analyser_20260210_143022.log
```

---

## Limitations

- Designed primarily for bacterial genomes
- GenBank files must contain CDS annotations for full functionality
- Memory usage scales with genome size

---

## What I Learned

- Structuring a Python project across multiple modules
- Using BioPython to parse and work with GenBank files
- Pandas for data manipulation and CSV export
- Matplotlib and Seaborn for scientific visualisation
- Argument parsing with `argparse`
- Logging across modules using Python's `logging` library
- Error handling with informative user-facing messages

---

## Learning Journey

This is my first bioinformatics project, built as a learning exercise with AI assistance (Claude AI). It represents my exploration of Python, bioinformatics, software architecture, and professional development practices. The code is functional and demonstrates what I've learned, though I will continue to improve it as my skills grow.


---

## Genome Data

Genome files are not included in this repository due to file size. The default genome used during development was the *E. coli* K-12 MG1655 complete genome, available from NCBI:

[https://www.ncbi.nlm.nih.gov/nuccore/U00096](https://www.ncbi.nlm.nih.gov/nuccore/U00096)

Download in GenBank format and place in `data/raw/`.

---

## Future Improvements

- [ ] **Progress bars** — Add `tqdm` progress bars to long-running phases like gene extraction and visualisation
- [ ] **Compressed file support** — Accept `.gz` genome files directly, since NCBI commonly distributes them in this format
- [ ] **Unit tests** — Write a pytest test suite
- [ ] **Config file support** — Allow analysis parameters (window size, step size, size bins) to be set via a YAML or TOML file instead of editing `config.py` directly
- [ ] **FASTA + GFF3 input** — Support the FASTA + GFF3 format pair, which many public datasets use instead of GenBank
- [ ] **Multi-genome batch mode** — Accept a directory of genome files and produce a comparative summary across all of them
- [ ] **Strand-separated density plots** — Extend the existing `split_genes_by_strand()` function to plot forward and reverse strand gene density separately
- [ ] **HTML report output** — Generate a self-contained HTML report with embedded interactive Plotly plots instead of static PNGs
- [ ] **Circular genome map** — Visualise genes, GC content, and strand distribution around a circular chromosome map

---

## Author

Jared Cambridge - January 2026

