# Copyright 2025, Battelle Energy Alliance, LLC, ALL RIGHTS RESERVED


import os
import argparse


def parse_CLI_arguments():
    """
    Parse command line arguments for filename input.
    
    Returns:
        str: The filename provided as a command line argument
    """
    parser = argparse.ArgumentParser(
        description="Process HydroBoost input files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python script.py Model_A
  python script.py Model_B
  python script.py Project_2
        """
    )
    
    parser.add_argument(
        "filename",
        type=str,
        help="Base filename (without path or extension)"
    )
    
    args = parser.parse_args()
    return args.filename


def generate_path(segments):
    """
    Create (if needed) and return the full path under:
      <repo_root>/Core_Models/HydroBoost/generated_data/<segment1>/<segment2>/...
    """
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    repo_root = os.path.abspath(os.path.join(script_dir, os.pardir, os.pardir))
    
    base = os.path.join(repo_root, "Core_Models", "HydroBoost", "generated_data")

    for seg in segments:
        base = os.path.join(base, seg)
        os.makedirs(base, exist_ok=True)

    return base

