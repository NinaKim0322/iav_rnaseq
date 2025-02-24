# RNA-seq Visualization Pipeline (Python)

This Python package provides a comprehensive pipeline for visualizing RNA-seq differential expression results from DESeq2, addressing the visualization errors you experienced with the R scripts.

## Features

- **Volcano plots**: Visualize gene significance vs. fold change for each comparison
- **MA plots**: Display mean expression vs. fold change with proper labeling
- **PCA plots**: Analyze sample clustering with confidence ellipses and statistical testing
- **Venn diagrams**: Compare differentially expressed gene sets between conditions
- **Heatmaps**: Visualize expression patterns of top differentially expressed genes
- **HTML Summary Report**: Generate a comprehensive report of all analysis results

## Installation

1. Clone this repository:
   ```
   git clone https://github.com/your-username/rna-seq-visualizer.git
   cd rna-seq-visualizer
   ```

2. Install the required dependencies:
   ```
   pip install -r requirements.txt
   ```

## Usage

1. Ensure your data files are in the current working directory:
   - DESeq2 output files:
     - `deseq2_male_3day.csv`
     - `deseq2_male_8day.csv`
     - `deseq2_female_3day.csv`
     - `deseq2_female_8day.csv`
   - `counts.csv` - Raw count matrix
   - `design.csv` - Sample metadata

2. Run the main script:
   ```
   python rna_seq_visualizer.py
   ```

3. Check the `results` directory for all generated visualizations and the summary report.

## Customization

You can customize the visualization parameters by modifying the function arguments in the main script. Here are some examples:

- Adjust volcano plot thresholds:
```python
process_volcano_plots(deseq2_files, fc_threshold=1.5, pval_threshold=0.01)
```

- Change MA plot y-axis limits:
```python
y_limits = {
    'male_3day': (-6, 6),
    'male_8day': (-8, 8),
    'female_3day': (-4, 4),
    'female_8day': (-10, 10)
}
process_ma_plots(deseq2_files, y_limits=y_limits)
```

## Output

The pipeline generates the following outputs in the `results` directory:

- **Volcano plots**: `volcano_DESeq2_*.png`
- **MA plots**: `MA_DESeq2_*.png`
- **PCA plot**: `PCA_by_group.png`
- **Venn diagrams**: `venn_*.png`
- **Heatmap**: `top_DEGs_expression_heatmap.png`
- **Gene lists**: `venn_list/*.csv`
- **HTML report**: `report/DESeq2_RNA_seq_analysis_report.html`

## Troubleshooting

If you encounter any issues:

1. Ensure all required files are present in the working directory
2. Check that all dependencies are installed correctly
3. Look for error messages in the console output
4. Adjust memory settings if dealing with large datasets:
   ```
   export PYTHONMEM=16G
   python rna_seq_visualizer.py
   ```

## Dependencies

- numpy
- pandas
- matplotlib
- seaborn
- scipy
- scikit-learn
- matplotlib-venn
- adjustText

## License

This project is licensed under the MIT License - see the LICENSE file for details.
