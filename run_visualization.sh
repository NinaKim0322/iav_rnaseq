#!/bin/bash
# RNA-seq Visualization Pipeline Runner Script
# This script helps run the Python visualization pipeline

# Display help information
show_help() {
    echo "RNA-seq Visualization Pipeline Runner"
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -h, --help           Display this help message"
    echo "  -d, --dir DIR        Set the working directory containing input files"
    echo "  -o, --output DIR     Set output directory (default: results)"
    echo "  -v, --verbose        Enable verbose output"
    echo "  -m, --memory SIZE    Set memory allocation (e.g., 16G)"
    echo "  --install            Install required dependencies first"
    echo ""
}

# Default settings
WORK_DIR=""
OUTPUT_DIR="results"
VERBOSE=""
MEMORY=""

# Parse command line arguments
while [ "$1" != "" ]; do
    case $1 in
        -h | --help )
            show_help
            exit 0
            ;;
        -d | --dir )
            shift
            WORK_DIR=$1
            ;;
        -o | --output )
            shift
            OUTPUT_DIR=$1
            ;;
        -v | --verbose )
            VERBOSE="--verbose"
            ;;
        -m | --memory )
            shift
            MEMORY=$1
            ;;
        --install )
            INSTALL_DEPS=1
            ;;
        * )
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
    shift
done

# Install dependencies if requested
if [ "$INSTALL_DEPS" = "1" ]; then
    echo "Installing required dependencies..."
    pip install -r requirements.txt
    if [ $? -ne 0 ]; then
        echo "Error installing dependencies. Please check your Python installation."
        exit 1
    fi
    echo "Dependencies installed successfully."
fi

# Set Python memory allocation if specified
if [ "$MEMORY" != "" ]; then
    export PYTHONMEM=$MEMORY
    echo "Set Python memory allocation to $MEMORY"
fi

# Construct the command
CMD="python rna_seq_visualizer.py"

if [ "$WORK_DIR" != "" ]; then
    CMD="$CMD --dir \"$WORK_DIR\""
fi

if [ "$OUTPUT_DIR" != "" ]; then
    CMD="$CMD --output \"$OUTPUT_DIR\""
fi

if [ "$VERBOSE" != "" ]; then
    CMD="$CMD $VERBOSE"
fi

# Print the command being run
echo "Running command: $CMD"
echo "-------------------------------------------"

# Run the visualization script
eval $CMD

# Check exit status
STATUS=$?
if [ $STATUS -ne 0 ]; then
    echo "Visualization failed with exit code $STATUS"
    echo "Please check the error messages above or see TROUBLESHOOTING.md for help."
    exit $STATUS
fi

exit 0