"""
Experimental data loading module
"""

from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd


def load_experiment_data(filepath: Path) -> np.ndarray:
    """
    Load experimental data

    Parameters
    ----------
    filepath : Path
        Path to the experimental data CSV file

    Returns
    -------
    np.ndarray, shape (n_points, 4)
        Experimental data [power_W, width_um, depth_um, area_um2]
        Note: if the area column is missing from the CSV, it will be filled with 0
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Experimental data file not found: {filepath}")

    df = pd.read_csv(filepath)

    # Check required columns (power, width, depth are required; area is optional)
    required_cols = ["power_W", "width_um", "depth_um"]
    optional_cols = ["area_um2"]

    # Attempt column name mapping
    alt_mapping = {
        "power_W": ["power", "Power", "P"],
        "width_um": ["width", "Width", "w"],
        "depth_um": ["depth", "Depth", "d"],
        "area_um2": ["area", "Area", "a"],
    }

    for col, alts in alt_mapping.items():
        if col not in df.columns:
            for alt in alts:
                if alt in df.columns:
                    df = df.rename(columns={alt: col})
                    break

    # Check required columns
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        raise ValueError(f"Experimental data is missing required columns: {missing_cols}")

    # If area column is missing, fill with 0
    if "area_um2" not in df.columns:
        print("Warning: area_um2 column is missing from the experimental data; filling with 0")
        df["area_um2"] = 0.0

    return df[["power_W", "width_um", "depth_um", "area_um2"]].values


def prepare_multivariate_data(filepath: Path) -> np.ndarray:
    """
    Prepare multivariate experimental data

    Parameters
    ----------
    filepath : Path
        Path to the experimental data CSV file

    Returns
    -------
    np.ndarray, shape (n_points, 4)
        Experimental data [power, width, depth, area]
    """
    print(f"Reading experimental data (multivariate format): {filepath}")
    data = load_experiment_data(filepath)
    print(f"Experimental data shape: {data.shape}, number of power points: {len(data)}\n")
    return data


def split_experiment_data(
    data: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Split experimental data into inputs and outputs

    Parameters
    ----------
    data : np.ndarray, shape (n, 4)
        Complete experimental data

    Returns
    -------
    x : np.ndarray, shape (n, 1)
        Input (power)
    y : np.ndarray, shape (n, 3)
        Output (width, depth, area)
    """
    x = data[:, 0:1]  # power
    y = data[:, 1:4]  # width, depth, area
    return x, y


def get_power_points(data: np.ndarray) -> np.ndarray:
    """
    Get the list of power points in the experimental data

    Parameters
    ----------
    data : np.ndarray
        Experimental data

    Returns
    -------
    np.ndarray
        List of unique power points
    """
    return np.unique(data[:, 0])
