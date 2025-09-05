#

# HydroBoost Model

# Main Authors:
# Jonghwan Kwon; Argonne National Laboratory; kwonj@anl.gov
# Carlos Josue Lopez; Argonne National Laboratory; clopezsalgado@anl.gov
# Alberto Grimaldi; Argonne National Laboratory; agrimaldi@anl.gov

# Current version: 2.0
# Last update: 07.31.2025

#

import os

def generate_path(segments):
    """
    Create (if needed) and return the full path under:
      <repo_root>/Core_Models/HydroBoost/generated_data/<segment1>/<segment2>/...
    where repo_root is C:/Users/agrimaldi/a-leaf-dev_HydroBoost_AG_version.
    """
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    repo_root = os.path.abspath(os.path.join(script_dir, os.pardir, os.pardir))
    
    base = os.path.join(repo_root, "Core_Models", "HydroBoost", "generated_data")

    for seg in segments:
        base = os.path.join(base, seg)
        os.makedirs(base, exist_ok=True)

    return base