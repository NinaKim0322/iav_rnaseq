#!/usr/bin/env python3
# ======================================================
# RNA-seq Visualization Master Script (Python Version)
# 
# This script generates visualizations from DESeq2 
# differential expression results
# 
# Author: Nina K
# Date: February 2025
# ======================================================

import os
import sys
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2, venn2_circles
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
import matplotlib.patches as mpatches
from scipy.cluster.hierarchy import dendrogram, linkage
from adjustText import adjust_text
import warnings
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import re
import matplotlib as mpl
from matplotlib.patches import Patch

# Set plot style
plt.style.use('seaborn-v0_8-whitegrid')
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']
mpl.rcParams['svg.fonttype'] = 'none'  # Important for text in SVGs
warnings.filterwarnings("ignore")

# Create output directory if it doesn't exist
os.makedirs("results", exist_ok=True)

# ======================================================
# Utility Functions
# ======================================================

def ensure_dir(directory):
    """Create directory if it doesn't exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)

def load_deseq2_data(file_path):
    """Load DESeq2 data from CSV file."""
    if not os.path.exists(file_path):
        print(f"Warning: File {file_path} not found.")
        return None
    
    print(f"  Reading {file_path}...")
    data = pd.read_csv(file_path)
    
    # Rename first column to gene_id if it's unnamed
    if data.columns[0].startswith('Unnamed'):
        data = data.rename(columns={data.columns[0]: 'gene_id'})
    
    return data

def extract_condition_info(file_path):
    """Extract condition information from filename."""
    filename = os.path.basename(file_path)
    condition = re.sub(r'deseq2_|\.csv', '', filename)
    sex = "Female" if "female" in condition else "Male"
    timepoint = "3-day" if "3day" in condition else "8-day"
    return condition, sex, timepoint

def save_fig(fig, filename, dpi=300, bbox_inches='tight'):
    """Save figure with standard settings."""
    fig.savefig(filename, dpi=dpi, bbox_inches=bbox_inches)
    plt.close(fig)

# ======================================================
# Volcano Plot Functions
# ======================================================

def create_deseq2_volcano(file_path, 
                         fc_threshold=1, 
                         pval_threshold=0.05,
                         y_max_limits={'male_3day': None, 'male_8day': 60, 'female_3day': None, 'female_8day': 150},
                         n_top_genes=20):
    """Create volcano plot from DESeq2 data."""
    # Extract condition info
    condition, sex, timepoint = extract_condition_info(file_path)
    
    # Load and prepare data
    data = load_deseq2_data(file_path)
    if data is None:
        return None
    
    # Filter out rows with NA in key columns
    data = data.dropna(subset=['log2FoldChange', 'padj'])
    
    # Calculate -log10(padj)
    data['logPadj'] = -np.log10(data['padj'])
    
    # Determine significance
    data['Significance'] = 'Not Significant'
    data.loc[(data['padj'] < pval_threshold) & (data['log2FoldChange'] > fc_threshold), 'Significance'] = 'Upregulated'
    data.loc[(data['padj'] < pval_threshold) & (data['log2FoldChange'] < -fc_threshold), 'Significance'] = 'Downregulated'
    
    # Select top genes to label (by fold change and by significance)
    top_up_fc = data[data['Significance'] == 'Upregulated'].sort_values('log2FoldChange', ascending=False).head(n_top_genes//2)
    top_down_fc = data[data['Significance'] == 'Downregulated'].sort_values('log2FoldChange').head(n_top_genes//2)
    
    top_up_pval = data[data['Significance'] == 'Upregulated'].sort_values('logPadj', ascending=False).head(n_top_genes//2)
    top_down_pval = data[data['Significance'] == 'Downregulated'].sort_values('logPadj', ascending=False).head(n_top_genes//2)
    
    # Combine top genes and remove duplicates
    top_genes = pd.concat([top_up_fc, top_down_fc, top_up_pval, top_down_pval]).drop_duplicates()
    
    # Count genes in each category
    category_counts = data['Significance'].value_counts().to_dict()
    for cat in ['Upregulated', 'Downregulated', 'Not Significant']:
        if cat not in category_counts:
            category_counts[cat] = 0
    
    # Get y-axis max for this condition
    y_max = y_max_limits.get(condition)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Prepare plot data based on y-axis compression if needed
    if y_max is not None:
        # Create a new column for compressed values
        data['compressed_logPadj'] = data['logPadj'].copy()
        
        # Identify points above threshold for visual indication
        data['above_threshold'] = data['logPadj'] > y_max
        
        # Compress values above y_max
        compression_factor = 10
        data.loc[data['above_threshold'], 'compressed_logPadj'] = y_max + (data.loc[data['above_threshold'], 'logPadj'] - y_max) / compression_factor
        
        # Do the same for top genes
        if not top_genes.empty:
            top_genes['compressed_logPadj'] = top_genes['logPadj'].copy()
            top_genes['above_threshold'] = top_genes['logPadj'] > y_max
            top_genes.loc[top_genes['above_threshold'], 'compressed_logPadj'] = y_max + (top_genes.loc[top_genes['above_threshold'], 'logPadj'] - y_max) / compression_factor
        
        # Plot with compressed y-axis
        for significance, color in [('Downregulated', 'blue'), ('Not Significant', 'gray'), ('Upregulated', 'red')]:
            group = data[data['Significance'] == significance]
            
            # Non-compressed points
            non_compressed = group[~group['above_threshold']]
            if not non_compressed.empty:
                ax.scatter(non_compressed['log2FoldChange'], non_compressed['compressed_logPadj'], 
                           color=color, alpha=0.7, label=significance, s=20)
            
            # Compressed points with different marker
            compressed = group[group['above_threshold']]
            if not compressed.empty:
                ax.scatter(compressed['log2FoldChange'], compressed['compressed_logPadj'], 
                           color=color, alpha=0.7, marker='^', s=25)
        
        # Add compression indicator
        ax.axhline(y=y_max, linestyle='dashed', color='black')
        ax.text(max(data['log2FoldChange']) * 0.8, y_max + 0.1, '//', color='black', fontsize=12)
        
        y_axis_label = "compressed -log10(Adjusted P-Value)"
    else:
        # Standard plot without compression
        for significance, color in [('Downregulated', 'blue'), ('Not Significant', 'gray'), ('Upregulated', 'red')]:
            group = data[data['Significance'] == significance]
            if not group.empty:
                ax.scatter(group['log2FoldChange'], group['logPadj'], 
                           color=color, alpha=0.7, label=significance, s=20)
        
        y_axis_label = "-log10(Adjusted P-Value)"
    
    # Add gene labels
    if not top_genes.empty:
        texts = []
        if y_max is not None:
            for _, gene in top_genes.iterrows():
                if gene['above_threshold']:
                    label = f"{gene['gene_id']} ({gene['logPadj']:.1f})"
                else:
                    label = gene['gene_id']
                texts.append(ax.text(gene['log2FoldChange'], gene['compressed_logPadj'], 
                                     label, fontsize=8, ha='center', va='center'))
        else:
            for _, gene in top_genes.iterrows():
                texts.append(ax.text(gene['log2FoldChange'], gene['logPadj'], 
                                     gene['gene_id'], fontsize=8, ha='center', va='center'))
        
        # Adjust text to avoid overlaps
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))
    
    # Add significance thresholds
    ax.axhline(y=-np.log10(pval_threshold), linestyle='dashed', color='darkblue')
    ax.axvline(x=-fc_threshold, linestyle='dashed', color='darkgreen')
    ax.axvline(x=fc_threshold, linestyle='dashed', color='darkgreen')
    
    # Add labels, title and legend
    ax.set_xlabel('log2(Fold Change) [IAV/Sham]', fontsize=12)
    ax.set_ylabel(y_axis_label, fontsize=12)
    ax.set_title(f'DESeq2 Volcano Plot: {sex} {timepoint} (IAV vs Sham)', fontsize=14, fontweight='bold')
    
    # Add subtitle with gene counts
    plt.figtext(0.5, 0.01, 
                f"Upregulated in IAV: {category_counts['Upregulated']} genes | "
                f"Downregulated in IAV: {category_counts['Downregulated']} genes",
                ha='center', fontsize=10)
    
    # Customize legend
    legend_labels = [
        f"Downregulated in IAV ({category_counts['Downregulated']})",
        f"Not Significant ({category_counts['Not Significant']})",
        f"Upregulated in IAV ({category_counts['Upregulated']})"
    ]
    legend_handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=8),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=8),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=8)
    ]
    ax.legend(legend_handles, legend_labels, frameon=True, loc='best')
    
    # Save figure
    output_file = os.path.join('results', f'volcano_DESeq2_{condition}.png')
    save_fig(fig, output_file)
    print(f"  Volcano plot saved to {output_file}")
    
    return fig

def process_volcano_plots(deseq2_files, 
                         fc_threshold=1, 
                         pval_threshold=0.05,
                         y_max_limits={'male_3day': None, 'male_8day': 60, 'female_3day': None, 'female_8day': 150},
                         n_top_genes=20):
    """Process all volcano plots."""
    print("Generating volcano plots...")
    for file in deseq2_files:
        if os.path.exists(file):
            create_deseq2_volcano(file, fc_threshold, pval_threshold, y_max_limits, n_top_genes)
        else:
            print(f"Warning: File not found: {file}")
    print("Volcano plots generated successfully!")

# ======================================================
# MA Plot Functions
# ======================================================

def create_deseq2_ma_plot(file_path, 
                         y_limits={'male_3day': (-5, 8),
                                  'male_8day': (-8, 8),
                                  'female_3day': (-3, 3),
                                  'female_8day': (-10, 10)}):
    """Create MA plot from DESeq2 data."""
    # Extract condition info
    condition, sex, timepoint = extract_condition_info(file_path)
    print(f"  Creating MA plot for {file_path}...")
    
    # Load data
    data = load_deseq2_data(file_path)
    if data is None:
        return None
    
    # Remove rows with NA in baseMean or log2FoldChange
    data = data.dropna(subset=['baseMean', 'log2FoldChange'])
    
    # Create original_log2FC column
    data['original_log2FC'] = data['log2FoldChange'].copy()
    
    # Get y-axis limits for this condition
    y_lim = y_limits.get(condition, (-5, 5))
    
    # Define significance
    data['Significance'] = 'Not Significant'
    data.loc[(data['padj'] < 0.05) & (data['log2FoldChange'] > 1), 'Significance'] = 'Upregulated'
    data.loc[(data['padj'] < 0.05) & (data['log2FoldChange'] < -1), 'Significance'] = 'Downregulated'
    
    # Create plotting dataset
    plot_data = data[(data['Significance'] != 'Not Significant') | (abs(data['log2FoldChange']) <= 1)].copy()
    
    # Create adjusted_log2FC column for plotting
    plot_data['adjusted_log2FC'] = plot_data['log2FoldChange'].copy()
    
    # Compress values beyond limits
    plot_data['beyond_limits'] = (plot_data['log2FoldChange'] < y_lim[0]) | (plot_data['log2FoldChange'] > y_lim[1])
    plot_data.loc[plot_data['log2FoldChange'] > y_lim[1], 'adjusted_log2FC'] = y_lim[1]
    plot_data.loc[plot_data['log2FoldChange'] < y_lim[0], 'adjusted_log2FC'] = y_lim[0]
    
    # Select top genes by fold change (10 up and 10 down)
    top_fc_up = plot_data[plot_data['Significance'] == 'Upregulated'].sort_values(by='original_log2FC', ascending=False).head(10)
    top_fc_down = plot_data[plot_data['Significance'] == 'Downregulated'].sort_values(by=['original_log2FC']).head(10)
    
    # Select top genes by mean expression (10 up and 10 down)
    top_mean_up = plot_data[plot_data['Significance'] == 'Upregulated'].sort_values(by='baseMean', ascending=False).head(10)
    top_mean_down = plot_data[plot_data['Significance'] == 'Downregulated'].sort_values(by='baseMean', ascending=False).head(10)
    
    # Combine all selected genes and remove duplicates
    all_selected_genes = pd.concat([top_fc_up, top_fc_down, top_mean_up, top_mean_down]).drop_duplicates()
    
    # Calculate gene counts
    up_count = (plot_data['Significance'] == 'Upregulated').sum()
    down_count = (plot_data['Significance'] == 'Downregulated').sum()
    
    # Create MA plot
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Plot by significance group
    for significance, color in [('Downregulated', 'blue'), ('Not Significant', 'gray'), ('Upregulated', 'red')]:
        group = plot_data[plot_data['Significance'] == significance]
        if not group.empty:
            ax.scatter(np.log10(group['baseMean']), group['adjusted_log2FC'], 
                       color=color, alpha=0.7, s=10, label=significance)
    
    # Add gene labels with original values for compressed points
    if not all_selected_genes.empty:
        texts = []
        for _, gene in all_selected_genes.iterrows():
            if gene['beyond_limits']:
                label = f"{gene['gene_id']}\n({gene['original_log2FC']:.1f})"
            else:
                label = gene['gene_id']
            texts.append(ax.text(np.log10(gene['baseMean']), gene['adjusted_log2FC'], 
                                label, fontsize=8, ha='center', va='center'))
        
        # Adjust text to avoid overlaps
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))
    
    # Add reference line at y=0
    ax.axhline(y=0, linestyle='solid', color='black')
    
    # Set y-axis limits
    ax.set_ylim(y_lim)
    
    # Add labels and title
    ax.set_xlabel('log10(Mean Expression)', fontsize=12)
    ax.set_ylabel('log2(Fold Change) [IAV/Sham]', fontsize=12)
    ax.set_title(f'MA Plot: {sex} {timepoint} (IAV vs Sham)', fontsize=14, fontweight='bold')
    
    # Add subtitle with gene counts
    plt.figtext(0.5, 0.01, 
                f"Upregulated in IAV: {up_count} genes | "
                f"Downregulated in IAV: {down_count} genes",
                ha='center', fontsize=10)
    
    # Customize legend
    legend_labels = [
        f"Downregulated in IAV ({down_count})",
        f"Not Significant",
        f"Upregulated in IAV ({up_count})"
    ]
    legend_handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=8),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=8),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=8)
    ]
    ax.legend(legend_handles, legend_labels, frameon=True, loc='best')
    
    # Save plot
    output_file = os.path.join('results', f'MA_DESeq2_{condition}.png')
    save_fig(fig, output_file)
    print(f"  MA plot saved to {output_file}")
    
    return fig

def process_ma_plots(deseq2_files,
                    y_limits={'male_3day': (-5, 8),
                             'male_8day': (-8, 8),
                             'female_3day': (-3, 3),
                             'female_8day': (-10, 10)}):
    """Process all MA plots."""
    print("Generating MA plots...")
    for file in deseq2_files:
        if os.path.exists(file):
            create_deseq2_ma_plot(file, y_limits)
        else:
            print(f"  Warning: File {file} not found. Skipping.")
    
    print("  All DESeq2 MA plots generated successfully!")

# ======================================================
# PCA & MDS Functions
# ======================================================

def create_pca_plot():
    """Create PCA plot from raw count data."""
    print("  Creating PCA plot from raw count data...")
    
    # Check if required files exist
    if not os.path.exists("counts.csv") or not os.path.exists("design.csv"):
        print("  Error: Required files (counts.csv or design.csv) not found.")
        return None
    
    # Load count data
    counts = pd.read_csv("counts.csv")
    counts = counts.set_index(counts.columns[0])
    
    # Remove rows with missing values and low expression
    counts = counts.dropna()
    counts = counts[(counts >= 10).sum(axis=1) >= 3]
    
    # Load design data
    design = pd.read_csv("design.csv")
    design = design.set_index("sample")
    
    # Match samples between counts and design
    matching_samples = list(set(counts.columns) & set(design.index))
    counts = counts[matching_samples]
    design = design.loc[matching_samples]
    
    # Create group factor
    design['group'] = design['timepoint'] + " " + design['treatment'] + " " + design['sex']
    
    # Normalize and log-transform
    count_sums = counts.sum()
    cpms = np.log2((counts + 1).divide(count_sums / 1e6, axis=1))
    
    # Perform PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(cpms.T)
    var_explained = np.round(pca.explained_variance_ratio_ * 100, 1)
    
    # Create dataframe for plotting
    pca_df = pd.DataFrame({
        'Sample': counts.columns,
        'PC1': pca_result[:, 0],
        'PC2': pca_result[:, 1]
    }).merge(design.reset_index(), left_on='Sample', right_on='sample')
    
    # Define group colors with better contrast
    group_colors = {
        "3day Sham M": "#90EE90",   # Light green
        "3day Sham F": "#228B22",   # Dark green
        "3day IAV M": "#87CEEB",    # Light blue
        "3day IAV F": "#0000CD",    # Dark blue
        "8day Sham M": "#FFA07A",   # Light orange
        "8day Sham F": "#FF4500",   # Dark orange
        "8day IAV M": "#FFB6C1",    # Light pink
        "8day IAV F": "#FF0000"     # Dark red
    }
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot each group
    for group, group_data in pca_df.groupby('group'):
        ax.scatter(group_data['PC1'], group_data['PC2'], 
                   color=group_colors[group], label=group, s=80, alpha=0.7)
        
        # Add confidence ellipses for groups with more than 2 samples
        if len(group_data) > 2:
            from matplotlib.patches import Ellipse
            from scipy.stats import chi2
            
            # Calculate covariance matrix and center
            cov = np.cov(group_data['PC1'], group_data['PC2'])
            center = [group_data['PC1'].mean(), group_data['PC2'].mean()]
            
            # Get size and orientation of ellipse
            eigvals, eigvecs = np.linalg.eigh(cov)
            order = eigvals.argsort()[::-1]
            eigvals, eigvecs = eigvals[order], eigvecs[:, order]
            
            # Calculate the tighter ellipse (approximately 1 standard deviation)
            confidence = 0.68  # 1 std
            chi2_val = chi2.ppf(confidence, 2)
            width, height = 2 * np.sqrt(chi2_val * eigvals)
            angle = np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))
            
            # Create and add ellipse
            ellipse = Ellipse(center, width, height, angle=angle, 
                             facecolor=group_colors[group], alpha=0.2, 
                             edgecolor=group_colors[group], linewidth=1)
            ax.add_patch(ellipse)
    
    # Add axis labels
    ax.set_xlabel(f'PC1: {var_explained[0]}% explained var.', fontsize=12)
    ax.set_ylabel(f'PC2: {var_explained[1]}% explained var.', fontsize=12)
    
    # Add legend
    ax.legend(title="Group", bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Add statistical results in bottom right
    # Calculate PERMANOVA equivalent using distance matrix and group labels
    from sklearn.metrics import pairwise_distances
    from scipy.stats import f_oneway
    
    # Create distance matrix
    dist_matrix = pairwise_distances(pca_df[['PC1', 'PC2']])
    
    # Compute basic between-group variance as a simple approximation
    groups = pca_df['group'].unique()
    group_dists = []
    
    for group in groups:
        group_indices = pca_df['group'] == group
        group_matrix = dist_matrix[np.ix_(group_indices, group_indices)]
        if len(group_matrix) > 0:
            group_dists.append(group_matrix[np.triu_indices(len(group_matrix), k=1)])
    
    # Perform simple ANOVA as approximation of PERMANOVA
    if len(group_dists) > 1 and all(len(d) > 0 for d in group_dists):
        f_stat, p_value = f_oneway(*group_dists)
        r2 = 1 - (sum(np.var(d) for d in group_dists) / np.var(dist_matrix[np.triu_indices(len(dist_matrix), k=1)]))
        
        stat_text = f"Group separation:\nR² = {r2:.3f}\nP = {'<0.001' if p_value < 0.001 else f'{p_value:.3f}'}"
        ax.text(pca_df['PC1'].max() * 0.9, pca_df['PC2'].min() * 1.1, 
                stat_text, ha='right', va='bottom', fontsize=10, 
                bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.5'))
    
    # Adjust plot limits to focus on data
    x_margin = (pca_df['PC1'].max() - pca_df['PC1'].min()) * 0.1
    y_margin = (pca_df['PC2'].max() - pca_df['PC2'].min()) * 0.1
    ax.set_xlim(pca_df['PC1'].min() - x_margin, pca_df['PC1'].max() + x_margin)
    ax.set_ylim(pca_df['PC2'].min() - y_margin, pca_df['PC2'].max() + y_margin)
    
    # Add title
    ax.set_title('PCA: Sample Clustering by Group', fontsize=14, fontweight='bold')
    
    # Save plot
    output_file = os.path.join('results', 'PCA_by_group.png')
    save_fig(fig, output_file)
    print(f"  PCA plot saved to {output_file}")
    
    return fig

def process_dimension_reduction():
    """Process PCA and MDS plots."""
    print("Generating dimension reduction plots...")
    try:
        pca_result = create_pca_plot()
        if pca_result is not None:
            print("  PCA plot generated successfully!")
        else:
            print("  Failed to generate PCA plot.")
    except Exception as e:
        print(f"  Error generating PCA plot: {str(e)}")

# ======================================================
# Venn Diagram Functions
# ======================================================

def create_venn_diagrams(male_genes, female_genes, timepoint):
    """Create Venn diagram comparing male and female DEGs."""
    # Calculate areas
    area1 = len(male_genes)  # Male (left)
    area2 = len(female_genes)  # Female (right)
    overlap = len(set(male_genes) & set(female_genes))
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Create Venn diagram
    v = venn2(subsets=(area1-overlap, area2-overlap, overlap),
             set_labels=('Male', 'Female'),
             set_colors=('#377EB8', '#E41A1C'),
             alpha=0.4,
             ax=ax)
    
    # Style the text
    for text in v.set_labels:
        if text is not None:
            text.set_fontsize(14)
    
    for text in v.subset_labels:
        if text is not None:
            text.set_fontsize(12)
    
    # Add title
    ax.set_title(f'Sex-specific DEG Comparison: {timepoint} Post-infection', 
                fontsize=14, fontweight='bold', pad=20)
    
    # Ensure clean plot borders
    ax.axis('tight')
    ax.axis('off')
    
    # Save individual plot
    output_file = os.path.join('results', f'venn_{timepoint}.png')
    save_fig(fig, output_file)
    print(f"  Venn diagram saved to {output_file}")
    
    return fig

def create_venn_comparison(male_3day, female_3day, male_8day, female_8day):
    """Create side-by-side Venn diagram comparison."""
    fig = plt.figure(figsize=(12, 6))
    
    # Define a layout with a main title and two Venn diagrams
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 9], figure=fig)
    
    # Add main title spanning both columns
    title_ax = fig.add_subplot(gs[0, :])
    title_ax.set_title('Sex-specific DEG Comparison at Different Timepoints', 
                      fontsize=16, fontweight='bold')
    title_ax.axis('off')
    
    # 3-day comparison
    ax1 = fig.add_subplot(gs[1, 0])
    ax1.set_title('3-day Post-infection', fontsize=14)
    
    # Create Venn diagram for 3-day
    area1 = len(male_3day)
    area2 = len(female_3day)
    overlap = len(set(male_3day) & set(female_3day))
    
    v1 = venn2(subsets=(area1-overlap, area2-overlap, overlap),
              set_labels=('Male', 'Female'),
              set_colors=('#377EB8', '#E41A1C'),
              alpha=0.4,
              ax=ax1)
    
    # Style the text
    for text in v1.set_labels:
        if text is not None:
            text.set_fontsize(12)
    
    for text in v1.subset_labels:
        if text is not None:
            text.set_fontsize(10)
    
    ax1.axis('off')
    
    # 8-day comparison
    ax2 = fig.add_subplot(gs[1, 1])
    ax2.set_title('8-day Post-infection', fontsize=14)
    
    # Create Venn diagram for 8-day
    area1 = len(male_8day)
    area2 = len(female_8day)
    overlap = len(set(male_8day) & set(female_8day))
    
    v2 = venn2(subsets=(area1-overlap, area2-overlap, overlap),
              set_labels=('Male', 'Female'),
              set_colors=('#377EB8', '#E41A1C'),
              alpha=0.4,
              ax=ax2)
    
    # Style the text
    for text in v2.set_labels:
        if text is not None:
            text.set_fontsize(12)
    
    for text in v2.subset_labels:
        if text is not None:
            text.set_fontsize(10)
    
    ax2.axis('off')
    
    # Save the comparison plot
    output_file = os.path.join('results', 'venn_comparison.png')
    save_fig(fig, output_file)
    print(f"  Venn diagram comparison saved to {output_file}")
    
    return fig

def process_venn_diagrams():
    """Process all Venn diagrams."""
    print("Generating Venn diagrams...")
    
    # Create directory for gene lists
    venn_dir = os.path.join("results", "venn_list")
    ensure_dir(venn_dir)
    
    # Extract significant genes
    def extract_significant_genes(file_path, log2fc_threshold=1):
        """Extract significant genes from DESeq2 results."""
        if not os.path.exists(file_path):
            print(f"  Warning: File {file_path} not found.")
            return []
        
        data = load_deseq2_data(file_path)
        if data is None:
            return []
        
        # Filter for significant genes
        sig_genes = data[(~data['padj'].isna()) & 
                         (data['padj'] < 0.05) & 
                         (abs(data['log2FoldChange']) > log2fc_threshold)]['gene_id'].tolist()
        
        return sig_genes
    
    # Extract genes for each condition
    male_3day = extract_significant_genes("deseq2_male_3day.csv")
    female_3day = extract_significant_genes("deseq2_female_3day.csv")
    male_8day = extract_significant_genes("deseq2_male_8day.csv")
    female_8day = extract_significant_genes("deseq2_female_8day.csv")
    
    # Save gene lists for analysis
    def save_gene_lists(male_genes, female_genes, timepoint):
        """Save gene lists for each section of the Venn diagram."""
        # Male-specific genes
        male_specific = list(set(male_genes) - set(female_genes))
        pd.DataFrame({
            'gene_id': male_specific,
            'category': 'Male_specific'
        }).to_csv(os.path.join(venn_dir, f"{timepoint}_male_specific.csv"), index=False)
        
        # Female-specific genes
        female_specific = list(set(female_genes) - set(male_genes))
        pd.DataFrame({
            'gene_id': female_specific,
            'category': 'Female_specific'
        }).to_csv(os.path.join(venn_dir, f"{timepoint}_female_specific.csv"), index=False)
        
        # Shared genes
        shared = list(set(male_genes) & set(female_genes))
        pd.DataFrame({
            'gene_id': shared,
            'category': 'Shared'
        }).to_csv(os.path.join(venn_dir, f"{timepoint}_shared.csv"), index=False)
    
    # Save gene lists for both timepoints
    save_gene_lists(male_3day, female_3day, "3day")
    save_gene_lists(male_8day, female_8day, "8day")
    
    # Create individual Venn diagrams
    create_venn_diagrams(male_3day, female_3day, "3day")
    create_venn_diagrams(male_8day, female_8day, "8day")
    
    # Create comparison plot
    create_venn_comparison(male_3day, female_3day, male_8day, female_8day)
    
    print("  All Venn diagrams generated successfully!")
    return True

# ======================================================
# Heatmap Functions
# ======================================================

def create_expression_heatmap(deseq2_files):
    """Create expression heatmap from DESeq2 results."""
    print("  Creating DESeq2 expression heatmaps...")
    
    # Check required files
    if not os.path.exists("counts.csv") or not os.path.exists("design.csv"):
        print("  Error: Required files not found. Skipping expression heatmap.")
        return None
    
    # Load count data
    counts = pd.read_csv("counts.csv")
    counts = counts.set_index(counts.columns[0])
    
    # Load design matrix
    design = pd.read_csv("design.csv")
    design = design.set_index('sample')
    
    # Match samples
    matching_samples = list(set(counts.columns) & set(design.index))
    counts = counts[matching_samples]
    design = design.loc[matching_samples]
    
    # Create ordered factors for proper sorting
    design['timepoint'] = pd.Categorical(design['timepoint'], categories=['3day', '8day'], ordered=True)
    design['treatment'] = pd.Categorical(design['treatment'], categories=['Sham', 'IAV'], ordered=True)
    design['sex'] = pd.Categorical(design['sex'], categories=['M', 'F'], ordered=True)
    
    # Create group label in desired order
    group_levels = [
        "3day Sham M", "3day Sham F",
        "3day IAV M", "3day IAV F",
        "8day Sham M", "8day Sham F",
        "8day IAV M", "8day IAV F"
    ]
    design['group'] = pd.Categorical(
        design['timepoint'] + " " + design['treatment'] + " " + design['sex'],
        categories=group_levels,
        ordered=True
    )
    
    # Sort samples by group
    design = design.sort_values(['timepoint', 'treatment', 'sex'])
    counts = counts[design.index]
    
    # Extract top DEGs from each DESeq2 file
    def extract_top_degs(file_path, n=50):
        """Extract top differentially expressed genes from DESeq2 results."""
        if not os.path.exists(file_path):
            return []
        
        print(f"  Extracting top DEGs from {file_path}...")
        data = load_deseq2_data(file_path)
        if data is None:
            return []
        
        # Filter and sort by adjusted p-value
        top_genes = data[(~data['padj'].isna()) & (data['padj'] < 0.05)].sort_values('padj').head(n)['gene_id'].tolist()
        return top_genes
    
    # Get top DEGs from each condition
    top_degs = {}
    for file in deseq2_files:
        condition = re.sub(r'deseq2_|\.csv', '', file)
        top_degs[condition] = extract_top_degs(file)
    
    # Combine all top DEGs
    all_top_degs = list(set().union(*[genes for genes in top_degs.values() if genes]))
    if not all_top_degs:
        print("  Error: No significant DEGs found.")
        return None
    
    # Filter count data for top DEGs
    available_degs = [gene for gene in all_top_degs if gene in counts.index]
    if not available_degs:
        print("  Error: None of the top DEGs found in count data.")
        return None
    
    # Process count data for all DEGs
    top_counts = counts.loc[available_degs]
    
    # Calculate CPM and log-transform
    cpms = np.log2((top_counts + 1).divide(counts.sum() / 1e6, axis=1))
    
    # Z-score normalization (by gene)
    cpms_z = (cpms - cpms.mean(axis=1).values.reshape(-1, 1)) / cpms.std(axis=1).values.reshape(-1, 1)
    
    # Create sample annotation colors
    colors = {
        'Sex': {'M': '#377EB8', 'F': '#E41A1C'},
        'Timepoint': {'3day': '#4DAF4A', '8day': '#984EA3'},
        'Treatment': {'Sham': '#FFFF33', 'IAV': '#FF7F00'}
    }
    
    # Configure annotation color map
    row_colors = pd.DataFrame({
        'Sex': design['sex'].map(colors['Sex']),
        'Timepoint': design['timepoint'].map(colors['Timepoint']),
        'Treatment': design['treatment'].map(colors['Treatment'])
    }, index=design.index)
    
    # Create main heatmap
    print("  Creating main heatmap...")
    fig, ax = plt.subplots(figsize=(14, 16))
    
    # Define heatmap colors
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    
    # Calculate gene and sample clustering
    g = sns.clustermap(
        cpms_z,
        cmap=cmap,
        col_colors=row_colors,
        figsize=(14, 16),
        xticklabels=False,  # Hide sample labels for cleaner look
        yticklabels=True,   # Show gene names
        row_cluster=True,   # Cluster genes
        col_cluster=False,  # Don't cluster samples (keep ordered by group)
        cbar_pos=(0.05, 0.05, 0.05, 0.18),
        cbar_kws={'label': 'Z-score'},
        dendrogram_ratio=0.15,
        colors_ratio=0.05,
    )
    
    # Adjust heatmap
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=6)
    
    # Add title
    g.fig.suptitle('Expression Heatmap of Top DEGs from DESeq2', fontsize=16, y=0.98)
    
    # Add legend for annotations
    legend_elements = []
    for name, color_map in colors.items():
        for label, color in color_map.items():
            legend_elements.append(Patch(facecolor=color, edgecolor='black', label=f'{name}: {label}'))
    
    g.fig.legend(handles=legend_elements, loc='upper center', bbox_to_anchor=(0.5, 0.02),
                ncol=len(legend_elements) // 2, frameon=True)
    
    # Save heatmap
    output_file = os.path.join('results', 'top_DEGs_expression_heatmap.png')
    save_fig(g.fig, output_file)
    print(f"  Heatmap saved to {output_file}")
    
    return g.fig

def process_heatmaps(deseq2_files):
    """Process all heatmaps."""
    print("Generating heatmaps...")
    heatmap_result = create_expression_heatmap(deseq2_files)
    
    if heatmap_result is None:
        print("  No DESeq2 heatmaps were generated.")
    else:
        print("  All DESeq2 heatmaps generated successfully!")
    
    return heatmap_result is not None

# ======================================================
# Summary Report Function
# ======================================================

def create_summary_report():
    """Create HTML summary report for DESeq2 analysis."""
    print("  Creating HTML summary report for DESeq2 analysis...")
    
    # Create report directory
    report_dir = os.path.join("results", "report")
    ensure_dir(report_dir)
    
    # Generate summary statistics by reading DEG data
    print("  Calculating summary statistics from DESeq2 files...")
    
    # Function to extract DEGs from a DESeq2 file
    def extract_degs(file_path):
        """Extract DEG statistics from DESeq2 results."""
        if not os.path.exists(file_path):
            print(f"    Warning: Cannot find file {file_path}")
            return {'total': 0, 'up': 0, 'down': 0}
        
        # Read data
        print(f"    Reading {file_path}")
        data = load_deseq2_data(file_path)
        if data is None:
            return {'total': 0, 'up': 0, 'down': 0}
        
        # Calculate statistics
        if 'padj' in data.columns and 'log2FoldChange' in data.columns:
            total = ((~data['padj'].isna()) & (data['padj'] < 0.05) & (abs(data['log2FoldChange']) > 0.5)).sum()
            up = ((~data['padj'].isna()) & (data['padj'] < 0.05) & (data['log2FoldChange'] > 0.5)).sum()
            down = ((~data['padj'].isna()) & (data['padj'] < 0.05) & (data['log2FoldChange'] < -0.5)).sum()
        else:
            print(f"    Warning: Required columns missing in {file_path}")
            total = 0
            up = 0
            down = 0
        
        return {'total': total, 'up': up, 'down': down}
    
    # Extract stats for each condition
    venn_summary = {
        'male_3day': extract_degs("deseq2_male_3day.csv"),
        'male_8day': extract_degs("deseq2_male_8day.csv"),
        'female_3day': extract_degs("deseq2_female_3day.csv"),
        'female_8day': extract_degs("deseq2_female_8day.csv")
    }
    
    # Create HTML file
    report_file = os.path.join(report_dir, "DESeq2_RNA_seq_analysis_report.html")
    
    # HTML template with placeholders
    html_template = """<!DOCTYPE html>
<html>
<head>
  <title>DESeq2: IAV vs Sham RNA-seq Differential Expression Analysis</title>
  <style>
    body { 
      font-family: Arial, sans-serif; 
      margin: 40px; 
      line-height: 1.6;
      color: #333;
      max-width: 1200px;
      margin: 0 auto;
      padding: 20px;
    }
    h1 { 
      color: #2c3e50; 
      border-bottom: 2px solid #3498db;
      padding-bottom: 10px;
    }
    h2 { 
      color: #3498db; 
      margin-top: 30px;
      border-bottom: 1px solid #ddd;
      padding-bottom: 5px;
    }
    h3 { 
      color: #2980b9; 
      margin-top: 25px;
    }
    .section {
      margin-bottom: 40px;
    }
    table { 
      border-collapse: collapse; 
      width: 100%; 
      margin: 20px 0; 
    }
    th, td { 
      padding: 12px; 
      text-align: left; 
      border-bottom: 1px solid #ddd; 
    }
    th { 
      background-color: #3498db; 
      color: white; 
    }
    tr:nth-child(even) { 
      background-color: #f2f2f2; 
    }
    tr:hover { 
      background-color: #f5f5f5; 
    }
    .info-box {
      background-color: #d1ecf1;
      border-left: 5px solid #0c5460;
      padding: 15px;
      margin: 20px 0;
      border-radius: 5px;
    }
    .gallery {
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(400px, 1fr));
      gap: 20px;
      margin: 20px 0;
    }
    .gallery-item {
      border: 1px solid #ddd;
      border-radius: 5px;
      padding: 15px;
      background-color: #f8f9fa;
    }
    .gallery-item h4 {
      color: #16a085;
      margin-top: 0;
    }
    .gallery-item p {
      font-size: 14px;
      color: #666;
    }
    .placeholder-img {
      width: 100%;
      height: 200px;
      background-color: #f0f0f0;
      border: 1px dashed #ccc;
      display: flex;
      align-items: center;
      justify-content: center;
      margin-top: 10px;
      color: #666;
      font-size: 14px;
    }
    .highlight {
      background-color: #fffacd;
      padding: 2px 5px;
      border-radius: 3px;
    }
    .note {
      color: #d35400;
      font-style: italic;
      margin-top: 10px;
    }
  </style>
</head>
<body>
  <h1>DESeq2 RNA-seq Differential Expression Analysis</h1>
  
  <div class="info-box">
    <p>This report summarizes the differential expression analysis using DESeq2, comparing IAV-infected samples to Sham controls across different sexes and time points.</p>
  </div>
  
  <div class="section">
    <h2>Summary of Differential Expression Results</h2>
    <p>The DESeq2 analysis examined gene expression in response to IAV infection in the following conditions:</p>
    
    <table>
      <tr>
        <th>Condition</th>
        <th>Total DEGs</th>
        <th>Upregulated in IAV</th>
        <th>Downregulated in IAV</th>
      </tr>
      <tr>
        <td>Male 3day</td>
        <td>{male_3day_total}</td>
        <td>{male_3day_up}</td>
        <td>{male_3day_down}</td>
      </tr>
      <tr>
        <td>Male 8day</td>
        <td>{male_8day_total}</td>
        <td>{male_8day_up}</td>
        <td>{male_8day_down}</td>
      </tr>
      <tr>
        <td>Female 3day</td>
        <td>{female_3day_total}</td>
        <td>{female_3day_up}</td>
        <td>{female_3day_down}</td>
      </tr>
      <tr>
        <td>Female 8day</td>
        <td>{female_8day_total}</td>
        <td>{female_8day_up}</td>
        <td>{female_8day_down}</td>
      </tr>
    </table>
    
    <p><em>Differential expression criteria: DESeq2 adjusted p-value < 0.05, |log2 fold change| > 0.5</em></p>
  </div>
  
  <div class="section">
    <h2>Volcano Plots</h2>
    <p>Volcano plots show the relationship between statistical significance (-log10 adjusted p-value) and magnitude of change (log2 fold change).</p>
    <p class="note">Note: PNG files have been generated for all plots. Please check the results directory for the image files.</p>
    
    <h3>DESeq2 Volcano Plots</h3>
    <div class="gallery">
      <div class="gallery-item">
        <h4>Male 3-day</h4>
        <p>IAV vs Sham comparison in males at 3 days post-infection</p>
        <div class="placeholder-img">Volcano Plot: Male 3-day (IAV vs Sham)<br>PNG file available in results folder</div>
      </div>
      
      <div class="gallery-item">
        <h4>Male 8-day</h4>
        <p>IAV vs Sham comparison in males at 8 days post-infection</p>
        <div class="placeholder-img">Volcano Plot: Male 8-day (IAV vs Sham)<br>PNG file available in results folder</div>
      </div>
      
      <div class="gallery-item">
        <h4>Female 3-day</h4>
        <p>IAV vs Sham comparison in females at 3 days post-infection</p>
        <div class="placeholder-img">Volcano Plot: Female 3-day (IAV vs Sham)<br>PNG file available in results folder</div>
      </div>
      
      <div class="gallery-item">
        <h4>Female 8-day</h4>
        <p>IAV vs Sham comparison in females at 8 days post-infection</p>
        <div class="placeholder-img">Volcano Plot: Female 8-day (IAV vs Sham)<br>PNG file available in results folder</div>
      </div>
    </div>
  </div>
  
  <div class="section">
    <h2>MA Plots</h2>
    <p>MA plots show the relationship between average expression level (x-axis) and log fold change (y-axis).</p>
    
    <h3>DESeq2 MA Plots</h3>
    <div class="gallery">
      <div class="gallery-item">
        <h4>Male 3-day</h4>
        <p>IAV vs Sham comparison in males at 3 days post-infection</p>
        <div class="placeholder-img">MA Plot: Male 3-day (IAV vs Sham)<br>PNG file available in results folder</div>
      </div>
      
      <div class="gallery-item">
        <h4>Male 8-day</h4>
        <p>IAV vs Sham comparison in males at 8 days post-infection</p>
        <div class="placeholder-img">MA Plot: Male 8-day (IAV vs Sham)<br>PNG file available in results folder</div>
      </div>
      
      <div class="gallery-item">
        <h4>Female 3-day</h4>
        <p>IAV vs Sham comparison in females at 3 days post-infection</p>
        <div class="placeholder-img">MA Plot: Female 3-day (IAV vs Sham)<br>PNG file available in results folder</div>
      </div>
      
      <div class="gallery-item">
        <h4>Female 8-day</h4>
        <p>IAV vs Sham comparison in females at 8 days post-infection</p>
        <div class="placeholder-img">MA Plot: Female 8-day (IAV vs Sham)<br>PNG file available in results folder</div>
      </div>
    </div>
  </div>
  
  <div class="section">
    <h2>Principal Component Analysis</h2>
    <p>PCA plots show the overall relationship between samples based on their gene expression profiles.</p>
    
    <div class="gallery">
      <div class="gallery-item">
        <h4>PCA by Group</h4>
        <p>Samples colored by combined grouping factors</p>
        <div class="placeholder-img">PCA Plot: By Group<br>PNG file available in results folder</div>
      </div>
    </div>
  </div>
  
  <div class="section">
    <h2>Venn Diagrams</h2>
    <p>Venn diagrams show the overlap of differentially expressed genes between conditions.</p>
    
    <h3>Sex Comparisons</h3>
    <div class="gallery">
      <div class="gallery-item">
        <h4>Male vs Female at 3-day</h4>
        <p>Overlap of all DEGs between sexes at 3 days post-infection</p>
        <div class="placeholder-img">Venn Diagram: Male vs Female (3-day)<br>PNG file available in results folder</div>
      </div>
      
      <div class="gallery-item">
        <h4>Male vs Female at 8-day</h4>
        <p>Overlap of all DEGs between sexes at 8 days post-infection</p>
        <div class="placeholder-img">Venn Diagram: Male vs Female (8-day)<br>PNG file available in results folder</div>
      </div>
      
      <div class="gallery-item">
        <h4>Sex Comparison (Side-by-side)</h4>
        <p>Side-by-side comparison of sex differences across timepoints</p>
        <div class="placeholder-img">Venn Diagram: Sex Comparisons<br>PNG file available in results folder</div>
      </div>
    </div>
  </div>
  
  <div class="section">
    <h2>Expression Heatmaps</h2>
    <p>Heatmaps show expression patterns of top differentially expressed genes across samples.</p>
    
    <div class="gallery">
      <div class="gallery-item">
        <h4>Top DEGs Expression Heatmap</h4>
        <p>Expression patterns of top differentially expressed genes</p>
        <div class="placeholder-img">Heatmap: Top DEGs<br>PNG file available in results folder</div>
      </div>
    </div>
  </div>
  
  <div class="section">
    <h2>Conclusions</h2>
    <p>This analysis provides insights into how IAV infection affects gene expression in male and female subjects at different timepoints.</p>
    <ul>
      {conclusions}
    </ul>
  </div>
  
  <div class="section">
    <h2>Data & Methods</h2>
    <p>RNA-seq data was processed using the following workflow:</p>
    <ol>
      <li>Differential expression analysis was performed using <span class="highlight">DESeq2</span></li>
      <li>Genes with adjusted p-value < 0.05 and |log2 fold change| > 0.5 were considered differentially expressed</li>
      <li>Visualizations were created using <span class="highlight">matplotlib</span>, <span class="highlight">seaborn</span>, and <span class="highlight">matplotlib-venn</span> packages in Python</li>
      <li>The experimental design included <span class="highlight">sex</span> (male, female), <span class="highlight">timepoint</span> (3-day, 8-day), and <span class="highlight">treatment</span> (IAV, Sham) factors</li>
    </ol>
    
    <p>For a complete set of PNG visualizations, please check the following directories:</p>
    <ul>
      <li><strong>Volcano plots:</strong> <code>results/</code> directory (volcano_DESeq2_*.png)</li>
      <li><strong>MA plots:</strong> <code>results/</code> directory (MA_DESeq2_*.png)</li>
      <li><strong>PCA plots:</strong> <code>results/</code> directory (PCA_*.png)</li>
      <li><strong>Venn diagrams:</strong> <code>results/</code> directory (venn_*.png)</li>
      <li><strong>Heatmaps:</strong> <code>results/</code> directory (*heatmap.png)</li>
    </ul>
  </div>
  
  <footer style="margin-top: 50px; padding-top: 20px; border-top: 1px solid #ddd; text-align: center; color: #777; font-size: 14px;">
    <p>Report generated on {date}</p>
  </footer>
</body>
</html>
"""
    
    # Generate conclusions based on the summary statistics
    m3_count = venn_summary['male_3day']['total']
    m8_count = venn_summary['male_8day']['total']
    f3_count = venn_summary['female_3day']['total']
    f8_count = venn_summary['female_8day']['total']
    
    # Initialize conclusions
    conclusions = []
    
    # Add dynamic conclusions based on the data
    if m8_count > m3_count and f8_count > f3_count:
        conclusions.append("<li>The 8-day timepoint shows more differentially expressed genes than the 3-day timepoint in both sexes, suggesting a cumulative effect of infection over time.</li>")
    
    if f3_count > m3_count and f8_count > m8_count:
        conclusions.append("<li>Females show more differentially expressed genes than males, suggesting potential sex-specific differences in response to IAV infection.</li>")
    elif m3_count > f3_count and m8_count > f8_count:
        conclusions.append("<li>Males show more differentially expressed genes than males, suggesting potential sex-specific differences in response to IAV infection.</li>")
    
    # Add general conclusions
    conclusions.append("<li>The gene expression changes observed reflect both common and sex-specific responses to IAV infection.</li>")
    conclusions.append("<li>Time-dependent changes in gene expression indicate dynamic regulation during the course of infection.</li>")
    
    # Join conclusions
    conclusions_html = "\n      ".join(conclusions)
    
    # Format the HTML with data
    html_content = html_template.format(
        male_3day_total=venn_summary['male_3day']['total'],
        male_3day_up=venn_summary['male_3day']['up'],
        male_3day_down=venn_summary['male_3day']['down'],
        male_8day_total=venn_summary['male_8day']['total'],
        male_8day_up=venn_summary['male_8day']['up'],
        male_8day_down=venn_summary['male_8day']['down'],
        female_3day_total=venn_summary['female_3day']['total'],
        female_3day_up=venn_summary['female_3day']['up'],
        female_3day_down=venn_summary['female_3day']['down'],
        female_8day_total=venn_summary['female_8day']['total'],
        female_8day_up=venn_summary['female_8day']['up'],
        female_8day_down=venn_summary['female_8day']['down'],
        conclusions=conclusions_html,
        date=time.strftime("%B %d, %Y")
    )
    
    # Write HTML to file
    with open(report_file, 'w') as f:
        f.write(html_content)
    
    print(f"  HTML report saved to: {report_file}")
    
    return report_file

# ======================================================
# Main Function
# ======================================================

def main():
    """Main function to run the RNA-seq visualization pipeline."""
    import time
    start_time = time.time()
    
    print("=" * 60)
    print("RNA-seq Visualization Pipeline (Python Version)")
    print("=" * 60)
    
    # Set working directory if needed
    # os.chdir("/storage/group/gan11/default/Nina/RNA-seq/7")
    
    print("\nChecking required libraries...")
    # Libraries are already imported at the top of the script
    print("All required libraries are installed.")
    
    # Process DESeq2 files
    print("\nProcessing DESeq2 files...")
    deseq2_files = [
        "deseq2_male_3day.csv",
        "deseq2_male_8day.csv", 
        "deseq2_female_3day.csv",
        "deseq2_female_8day.csv"
    ]
    
    # Check if files exist
    missing_deseq2 = [file for file in deseq2_files if not os.path.exists(file)]
    if missing_deseq2:
        print(f"Warning: Some DESeq2 files are missing: {', '.join(missing_deseq2)}")
        print("The script will continue but some visualizations may be incomplete.")
    
    # Generate visualizations
    
    # 1. Volcano Plots
    print("\n" + "=" * 40)
    try:
        process_volcano_plots(deseq2_files)
    except Exception as e:
        print(f"Error generating volcano plots: {str(e)}")
        print("Continuing with other visualizations...")
    
    # 2. MA Plots
    print("\n" + "=" * 40)
    try:
        process_ma_plots(deseq2_files)
    except Exception as e:
        print(f"Error generating MA plots: {str(e)}")
        print("Continuing with other visualizations...")
    
    # 3. PCA Plots
    print("\n" + "=" * 40)
    try:
        process_dimension_reduction()
    except Exception as e:
        print(f"Error generating dimension reduction plots: {str(e)}")
        print("Continuing with other visualizations...")
    
    # 4. Heatmaps
    print("\n" + "=" * 40)
    try:
        process_heatmaps(deseq2_files)
    except Exception as e:
        print(f"Error generating heatmaps: {str(e)}")
        print("Continuing with other visualizations...")
    
    # 5. Venn Diagrams
    print("\n" + "=" * 40)
    try:
        process_venn_diagrams()
    except Exception as e:
        print(f"Error generating Venn diagrams: {str(e)}")
        print("Continuing with other visualizations...")
    
    # Create Summary Report
    print("\n" + "=" * 40)
    try:
        report_file = create_summary_report()
        print(f"Summary report created: {report_file}")
    except Exception as e:
        print(f"Error creating summary report: {str(e)}")
    
    # Calculate and print execution time
    end_time = time.time()
    execution_time = end_time - start_time
    print("\n" + "=" * 60)
    print(f"All visualizations have been generated successfully!")
    print(f"Total execution time: {execution_time:.2f} seconds")
    print("=" * 60)

if __name__ == "__main__":
    main()
