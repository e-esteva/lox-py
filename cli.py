import argparse
import sys
import os
from .core import process_all_pairs

def main(argv=None):
    """CLI entrypoint for loxcode-align."""
    parser = argparse.ArgumentParser(
        description="LoxCode Alignment Tool - Process paired-end sequencing LoxCode data"
    )
    parser.add_argument("input_folder", help="Input folder containing FASTQ files")
    parser.add_argument("r1_pattern", help="Pattern for R1 files (e.g., 'R1_001.fastq')")
    parser.add_argument("r2_pattern", help="Pattern for R2 files (e.g., 'R2_001.fastq')")
    parser.add_argument("output_folder", help="Output folder for results")

    args = parser.parse_args(argv)

    if not os.path.exists(args.output_folder):
        print(f"Output folder not found. Please create it first: mkdir {args.output_folder}")
        sys.exit(1)

    try:
        process_all_pairs(args.input_folder, args.r1_pattern, args.r2_pattern, args.output_folder)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

