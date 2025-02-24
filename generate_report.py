#!/usr/bin/env python3
# ======================================================
# RNA-seq HTML Report Generator
# 
# This script generates an HTML summary report for RNA-seq
# visualizations created by rna_seq_visualizer.py
# 
# Author: Nina K
# Date: February 2025
# ======================================================

import os
import sys
import argparse
import datetime
import pandas as pd
import re

def ensure_dir(directory):
    """Create directory if it doesn't exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)
        print(f"  Created directory: {directory}")

def load_deseq2_data(file_path):
    """Load DESeq2 data from CSV file with error handling."""
    if not os.path.exists(file_path):
        print(f"  Warning: File {file_path} not found.")
        return None
    
    print(f"  Reading {file_path}...")
    try:
        data = pd.read_csv(file_path)
        
        # Rename first column to gene_id if it's unnamed
        if data.columns[0].startswith('Unnamed'):
            data = data.rename(columns={data.columns[0]: 'gene_id'})
        
        # Check for required columns
        required_columns = ['gene_id', 'log2FoldChange', 'padj']
        missing_columns = [col for col in required_columns if col not in data.columns]
        
        if missing_columns:
            print(f"  Warning: Required columns missing in {file_path}: {', '.join(missing_columns)}")
            return None
            
        return data
        
    except Exception as e:
        print(f"  Error loading {file_path}: {str(e)}")
        return None

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

def create_summary_report(output_dir):
    """Create HTML summary report for DESeq2 analysis."""
    print("Generating HTML summary report for DESeq2 analysis...")
    
    # Create report directory
    report_dir = os.path.join(output_dir, "report")
    ensure_dir(report_dir)
    
    # Generate summary statistics by reading DEG data
    print("Calculating summary statistics from DESeq2 files...")
    
    # Extract stats for each condition
    venn_summary = {
        'male_3day': extract_degs("deseq2_male_3day.csv"),
        'male_8day': extract_degs("deseq2_male_8day.csv"),
        'female_3day': extract_degs("deseq2_female_3day.csv"),
        'female_8day': extract_degs("deseq2_female_8day.csv")
    }
    
    # Create HTML file
    report_file = os.path.join(report_dir, "DESeq2_RNA_seq_analysis_report.html")
    
    # Initialize conclusions
    conclusions = []
    
    # Add dynamic conclusions based on the data
    m3_count = venn_summary['male_3day']['total']
    m8_count = venn_summary['male_8day']['total']
    f3_count = venn_summary['female_3day']['total']
    f8_count = venn_summary['female_8day']['total']
    
    # Add dynamic conclusions based on the data
    if m8_count > m3_count and f8_count > f3_count:
        conclusions.append("<li>The 8-day timepoint shows more differentially expressed genes than the 3-day timepoint in both sexes, suggesting a cumulative effect of infection over time.</li>")
    
    if f3_count > m3_count and f8_count > m8_count:
        conclusions.append("<li>Females show more differentially expressed genes than males, suggesting potential sex-specific differences in response to IAV infection.</li>")
    elif m3_count > f3_count and m8_count > f8_count:
        conclusions.append("<li>Males show more differentially expressed genes than females, suggesting potential sex-specific differences in response to IAV infection.</li>")
    
    # Add general conclusions
    conclusions.append("<li>The gene expression changes observed reflect both common and sex-specific responses to IAV infection.</li>")
    conclusions.append("<li>Time-dependent changes in gene expression indicate dynamic regulation during the course of infection.</li>")
    
    # Join conclusions
    conclusions_html = "\n      ".join(conclusions)
    
    # Format the HTML with data
    html_content = f"""<!DOCTYPE html>
<html>
<head>
  <title>DESeq2: IAV vs Sham RNA-seq Differential Expression Analysis</title>
  <style>
    body {{ 
      font-family: Arial, sans-serif; 
      margin: 40px; 
      line-height: 1.6;
      color: #333;
      max-width: 1200px;
      margin: 0 auto;
      padding: 20px;
    }}
    h1 {{ 
      color: #2c3e50; 
      border-bottom: 2px solid #3498db;
      padding-bottom: 10px;
    }}
    h2 {{ 
      color: #3498db; 
      margin-top: 30px;
      border-bottom: 1px solid #ddd;
      padding-bottom: 5px;
    }}
    h3 {{ 
      color: #2980b9; 
      margin-top: 25px;
    }}
    .section {{
      margin-bottom: 40px;
    }}
    table {{ 
      border-collapse: collapse; 
      width: 100%; 
      margin: 20px 0; 
    }}
    th, td {{ 
      padding: 12px; 
      text-align: left; 
      border-bottom: 1px solid #ddd; 
    }}
    th {{ 
      background-color: #3498db; 
      color: white; 
    }}
    tr:nth-child(even) {{ 
      background-color: #f2f2f2; 
    }}
    tr:hover {{ 
      background-color: #f5f5f5; 
    }}
    .info-box {{
      background-color: #d1ecf1;
      border-left: 5px solid #0c5460;
      padding: 15px;
      margin: 20px 0;
      border-radius: 5px;
    }}
    .gallery {{
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(400px, 1fr));
      gap: 20px;
      margin: 20px 0;
    }}
    .gallery-item {{
      border: 1px solid #ddd;
      border-radius: 5px;
      padding: 15px;
      background-color: #f8f9fa;
    }}
    .gallery-item h4 {{
      color: #16a085;
      margin-top: 0;
    }}
    .gallery-item p {{
      font-size: 14px;
      color: #666;
    }}
    .placeholder-img {{
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
    }}
    .highlight {{
      background-color: #fffacd;
      padding: 2px 5px;
      border-radius: 3px;
    }}
    .note {{
      color: #d35400;
      font-style: italic;
      margin-top: 10px;
    }}
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
        <td>{venn_summary['male_3day']['total']}</td>
        <td>{venn_summary['male_3day']['up']}</td>
        <td>{venn_summary['male_3day']['down']}</td>
      </tr>
      <tr>
        <td>Male 8day</td>
        <td>{venn_summary['male_8day']['total']}</td>
        <td>{venn_summary['male_8day']['up']}</td>
        <td>{venn_summary['male_8day']['down']}</td>
      </tr>
      <tr>
        <td>Female 3day</td>
        <td>{venn_summary['female_3day']['total']}</td>
        <td>{venn_summary['female_3day']['up']}</td>
        <td>{venn_summary['female_3day']['down']}</td>
      </tr>
      <tr>
        <td>Female 8day</td>
        <td>{venn_summary['female_8day']['total']}</td>
        <td>{venn_summary['female_8day']['up']}</td>
        <td>{venn_summary['female_8day']['down']}</td>
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
      {conclusions_html}
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
      <li><strong>Volcano plots:</strong> <code>{output_dir}/</code> directory (volcano_DESeq2_*.png)</li>
      <li><strong>MA plots:</strong> <code>{output_dir}/</code> directory (MA_DESeq2_*.png)</li>
      <li><strong>PCA plots:</strong> <code>{output_dir}/</code> directory (PCA_*.png)</li>
      <li><strong>Venn diagrams:</strong> <code>{output_dir}/</code> directory (venn_*.png)</li>
      <li><strong>Heatmaps:</strong> <code>{output_dir}/</code> directory (*heatmap.png)</li>
    </ul>
  </div>
  
  <footer style="margin-top: 50px; padding-top: 20px; border-top: 1px solid #ddd; text-align: center; color: #777; font-size: 14px;">
    <p>Report generated on {datetime.datetime.now().strftime("%B %d, %Y")}</p>
  </footer>
</body>
</html>"""
    
    # Write HTML to file
    with open(report_file, 'w') as f:
        f.write(html_content)
    
    print(f"HTML report saved to: {report_file}")
    
    return report_file

def main():
    """Main function to run the HTML report generator."""
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='RNA-seq HTML Report Generator')
    parser.add_argument('--dir', type=str, default='.', help='Working directory containing input files')
    parser.add_argument('--output', type=str, default='results', help='Output directory for results')
    args = parser.parse_args()
    
    # Set the working directory if specified
    if args.dir != '.':
        try:
            os.chdir(args.dir)
            print(f"Changed working directory to: {args.dir}")
        except Exception as e:
            print(f"Error changing to directory {args.dir}: {str(e)}")
            print("Using current directory instead.")
    
    # Check for required files
    deseq2_files = [
        "deseq2_male_3day.csv",
        "deseq2_male_8day.csv", 
        "deseq2_female_3day.csv",
        "deseq2_female_8day.csv"
    ]
    
    missing_files = [file for file in deseq2_files if not os.path.exists(file)]
    if missing_files:
        print(f"Warning: Some required files are missing: {', '.join(missing_files)}")
        print("The report may be incomplete.")
    
    # Create output directory
    ensure_dir(args.output)
    
    # Generate the HTML report
    try:
        report_file = create_summary_report(args.output)
        print(f"\nReport generation completed successfully!")
        print(f"HTML report saved to: {report_file}")
    except Exception as e:
        print(f"Error generating report: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()