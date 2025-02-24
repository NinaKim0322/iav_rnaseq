# Remove the HTML report generation functionality from rna_seq_visualizer.py

def main():
    """Main function to run the RNA-seq visualization pipeline."""
    start_time = time.time()
    
    print("=" * 60)
    print("RNA-seq Visualization Pipeline (Python Version)")
    print("=" * 60)
    
    # Set working directory if needed
    # os.chdir("/storage/group/gan11/default/Nina/RNA-seq/7")
    
    # Check if we're running in verbose mode
    import argparse
    parser = argparse.ArgumentParser(description='RNA-seq Visualization Pipeline')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose output')
    parser.add_argument('--dir', type=str, help='Working directory containing input files')
    parser.add_argument('--output', type=str, default='results', help='Output directory for results')
    args = parser.parse_args()
    
    # Set the working directory if specified
    if args.dir:
        try:
            os.chdir(args.dir)
            print(f"\nChanged working directory to: {args.dir}")
        except Exception as e:
            print(f"\nError changing to directory {args.dir}: {str(e)}")
            print("Using current directory instead.")
    
    # Create output directory
    output_dir = args.output
    ensure_dir(output_dir)
    
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
    
    # Print input file information
    if args.verbose:
        print("\nInput files:")
        for file in deseq2_files:
            if os.path.exists(file):
                file_size = os.path.getsize(file) / 1024  # KB
                print(f"  {file}: {file_size:.1f} KB")
            else:
                print(f"  {file}: Not found")
        
        if os.path.exists("counts.csv"):
            count_size = os.path.getsize("counts.csv") / (1024 * 1024)  # MB
            print(f"  counts.csv: {count_size:.1f} MB")
        else:
            print("  counts.csv: Not found")
            
        if os.path.exists("design.csv"):
            design_size = os.path.getsize("design.csv") / 1024  # KB
            print(f"  design.csv: {design_size:.1f} KB")
        else:
            print("  design.csv: Not found")
    
    # Generate visualizations
    results_summary = {
        'volcano_plots': False,
        'ma_plots': False,
        'pca_plots': False,
        'heatmaps': False,
        'venn_diagrams': False
    }
    
    # 1. Volcano Plots
    print("\n" + "=" * 40)
    try:
        process_volcano_plots(deseq2_files)
        results_summary['volcano_plots'] = True
    except Exception as e:
        print(f"Error generating volcano plots: {str(e)}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        print("Continuing with other visualizations...")
    
    # 2. MA Plots
    print("\n" + "=" * 40)
    try:
        process_ma_plots(deseq2_files)
        results_summary['ma_plots'] = True
    except Exception as e:
        print(f"Error generating MA plots: {str(e)}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        print("Continuing with other visualizations...")
    
    # 3. PCA Plots
    print("\n" + "=" * 40)
    try:
        process_dimension_reduction()
        results_summary['pca_plots'] = True
    except Exception as e:
        print(f"Error generating dimension reduction plots: {str(e)}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        print("Continuing with other visualizations...")
    
    # 4. Heatmaps
    print("\n" + "=" * 40)
    try:
        process_heatmaps(deseq2_files)
        results_summary['heatmaps'] = True
    except Exception as e:
        print(f"Error generating heatmaps: {str(e)}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        print("Continuing with other visualizations...")
    
    # 5. Venn Diagrams
    print("\n" + "=" * 40)
    try:
        process_venn_diagrams()
        results_summary['venn_diagrams'] = True
    except Exception as e:
        print(f"Error generating Venn diagrams: {str(e)}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        print("Continuing with other visualizations...")
    
    # Calculate and print execution time
    end_time = time.time()
    execution_time = end_time - start_time
    
    # Print summary of results
    print("\n" + "=" * 60)
    print("RESULTS SUMMARY:")
    print("-" * 60)
    for key, success in results_summary.items():
        status = "✓ Success" if success else "✗ Failed"
        print(f"{key.replace('_', ' ').title():20}: {status}")
    print("-" * 60)
    
    # Print overall success message
    successful = sum(1 for success in results_summary.values() if success)
    total = len(results_summary)
    
    print(f"\nCompleted {successful}/{total} visualization tasks.")
    print(f"Total execution time: {execution_time:.2f} seconds")
    print(f"Results saved to: {os.path.abspath(output_dir)}")
    print("=" * 60)
    print("\nTo generate an HTML report, run: python generate_report.py")