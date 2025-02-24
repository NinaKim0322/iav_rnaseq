# Troubleshooting Guide for RNA-seq Visualization Pipeline

This guide helps you resolve common issues that might arise when running the RNA-seq visualization pipeline.

## Font Warnings

### Problem:
```
findfont: Generic family 'sans-serif' not found because none of the following families were found: Arial
```

### Solution:
These warnings have been suppressed in the latest version. If you're still seeing them:

1. The code now uses 'DejaVu Sans' by default, which should be available on most systems
2. You can modify the font settings in the script:
   ```python
   # Set a different font family
   mpl.rcParams['font.family'] = 'monospace'  # Try this alternative
   ```
3. Install additional fonts on your system:
   - **Ubuntu/Debian**: `sudo apt-get install fonts-dejavu`
   - **CentOS/RHEL**: `sudo yum install dejavu-sans-fonts`
   - **macOS**: Install with Homebrew: `brew cask install font-dejavu-sans`

## Heatmap Error

### Problem:
```
Error generating heatmaps: unsupported operand type(s) for +: 'Categorical' and 'str'
```

### Solution:
This error has been fixed by properly converting categorical columns to strings before concatenation:

```python
# Convert categorical columns to strings before concatenating
design['group'] = (design['timepoint'].astype(str) + " " + 
                  design['treatment'].astype(str) + " " + 
                  design['sex'].astype(str))
```

## Missing Required Package

### Problem:
```
ImportError: No module named matplotlib_venn
```

### Solution:
1. Install the missing package:
   ```
   pip install matplotlib-venn
   ```
2. Ensure all required packages are installed:
   ```
   pip install -r requirements.txt
   ```

## File Not Found Error

### Problem:
```
FileNotFoundError: [Errno 2] No such file or directory: 'counts.csv'
```

### Solution:
1. Check that all required data files are in the current working directory:
   - DESeq2 results files: `deseq2_male_3day.csv`, `deseq2_male_8day.csv`, etc.
   - Raw count data: `counts.csv`
   - Sample metadata: `design.csv`
   
2. Verify the file paths:
   ```python
   import os
   print(os.getcwd())  # Print current working directory
   print(os.listdir('.'))  # List files in current directory
   ```

3. Set the working directory explicitly if needed:
   ```python
   os.chdir("/path/to/your/data/directory")
   ```

## Memory Issues with Large Datasets

### Problem:
The script crashes or becomes extremely slow with large datasets.

### Solution:
1. Increase Python's memory allocation:
   ```bash
   export PYTHONMEM=16G  # Allocate 16GB of RAM
   ```

2. Process data in smaller chunks by modifying the script:
   ```python
   # Example: Limit the number of genes in heatmap
   def create_expression_heatmap(deseq2_files, max_genes=1000):
       # Limit the number of genes
       all_top_degs = all_top_degs[:max_genes]
   ```

## Plot Quality Issues

### Problem:
Plots appear blurry or with poor resolution.

### Solution:
Increase the DPI (dots per inch) setting in the save_fig function:

```python
def save_fig(fig, filename, dpi=600, bbox_inches='tight'):
    """Save figure with higher resolution"""
    fig.savefig(filename, dpi=dpi, bbox_inches=bbox_inches)
    plt.close(fig)
```

## Python Version Compatibility

### Problem:
Script fails with syntax errors or package incompatibilities.

### Solution:
This script was developed for Python 3.6+. If you're using an older version:

1. Verify your Python version:
   ```bash
   python --version
   ```

2. Install a compatible Python version or update your existing installation.

3. Consider creating a virtual environment with the correct Python version:
   ```bash
   # Create and activate a Python 3.8 virtual environment
   python3.8 -m venv env
   source env/bin/activate  # On Linux/Mac
   # or
   env\Scripts\activate  # On Windows
   ```
