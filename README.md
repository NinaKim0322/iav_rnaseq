# RNA-seq Visualization Pipeline (Python)

This Python package provides a comprehensive pipeline for visualizing RNA-seq differential expression results from DESeq2, addressing the visualization errors you experienced with the R scripts. The visualization and reporting functionality have been separated into two different scripts for better modularity.

## Features

- **Volcano plots**: Visualize gene significance vs. fold change for each comparison
- **MA plots**: Display mean expression vs. fold change with proper labeling
- **PCA plots**: Analyze sample clustering with confidence ellipses and statistical testing
- **Venn diagrams**: Compare differentially expressed gene sets between conditions
- **Heatmaps**: Visualize expression patterns of top differentially expressed genes
- **HTML Summary Report**: Generate a comprehensive report of all analysis results (separate script)

## File Structure

- `rna_seq_visualizer.py` - Main script for generating visualizations
- `generate_report.py` - Script for creating the HTML summary report
- `requirements.txt` - List of required Python packages
- `environment.yml` - Conda environment configuration
- `run_visualization.sh` - Bash script to simplify execution

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/your-username/rna-seq-visualizer.git
