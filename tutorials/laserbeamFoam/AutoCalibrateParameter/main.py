#!/usr/bin/env python3
"""
Tutorial entry point for AutoCalibrateParameter.

This script calibrates laser beam welding simulation parameters against
experimental meltpool measurements using Bayesian or gradient-based
optimisation.

Usage:
    python main.py                        # default: Bayesian with config.yaml
    python main.py --method gradient      # gradient-based optimisation
    python main.py --help                 # show all options
"""
import sys
from pathlib import Path

# Add applications/scripts to Python path so we can import AutoCalibrateParameter
# Structure: LaserbeamFoam/tutorials/laserbeamFoam/AutoCalibrateParameter/main.py
repo_root = Path(__file__).resolve().parents[3]  # -> LaserbeamFoam
scripts_path = repo_root / 'applications' / 'scripts'
if str(scripts_path) not in sys.path:
    sys.path.insert(0, str(scripts_path))

print("Loading AutoCalibrateParameter…", flush=True)
try:
    from AutoCalibrateParameter.run_cli import main
except ImportError as e:
    print(f"Error: Could not import AutoCalibrateParameter from {scripts_path}")
    print(f"Details: {e}")
    print()
    print("Please install required packages first:")
    print("  conda env create -f ../../../applications/scripts/postProcessing/environment.yml")
    print("  conda activate meltpool-postproc")
    print("  pip install scikit-optimize scipy pyyaml")
    sys.exit(1)

if __name__ == "__main__":
    sys.exit(main())
