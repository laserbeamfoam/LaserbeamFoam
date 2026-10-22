"""
Synthetic data management module
"""

from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

try:
    from pyDOE import lhs
except ImportError:
    lhs = None


def load_synthetic_data(filepath: Path) -> np.ndarray:
    """
    Load synthetic data

    Parameters
    ----------
    filepath : Path
        Path to the synthetic data file

    Returns
    -------
    np.ndarray, shape (n_samples, 9)
        [power, sigma, marangoni, substrate_temp, absorptivity, recoilCoeff, radius_flavour, width, depth, area]
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Synthetic data file not found: {filepath}")

    data = np.loadtxt(filepath)
    print(f"Loaded synthetic data: {filepath}, shape {data.shape}")
    return data


def save_synthetic_data(
    filepath: Path,
    data: np.ndarray,
    header: str = "power sigma marangoni substrate_temp absorptivity recoilCoeff radius_flavour width depth area",
) -> None:
    """
    Save synthetic data

    Parameters
    ----------
    filepath : Path
        Output file path
    data : np.ndarray
        Synthetic data
    header : str
        File header comment
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(filepath, data, header=header)
    print(f"Synthetic data saved to: {filepath}, shape: {data.shape}")


def append_synthetic_data(filepath: Path, new_data: np.ndarray) -> np.ndarray:
    """
    Append new data to the synthetic data file

    Parameters
    ----------
    filepath : Path
        Path to the synthetic data file
    new_data : np.ndarray
        New data to append

    Returns
    -------
    np.ndarray
        Updated complete dataset
    """
    import shutil

    existing_data = load_synthetic_data(filepath)
    updated_data = np.vstack([existing_data, new_data])

    # Back up the original file
    backup_path = filepath.with_suffix(filepath.suffix + ".backup")
    shutil.copy(filepath, backup_path)

    # Save the updated data
    save_synthetic_data(filepath, updated_data)
    print(f"Data updated: {filepath}, new shape {updated_data.shape}")

    return updated_data


def generate_lhs_samples(
    n_samples: int,
    param_bounds: List[Tuple[float, float]],
    power_range: Tuple[float, float] = (140.0, 260.0),
    seed: Optional[int] = 42,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate samples using Latin hypercube sampling

    Parameters
    ----------
    n_samples : int
        Number of samples
    param_bounds : list of tuple
        Parameter bounds list [(min, max), ...]
    power_range : tuple
        Power range (min, max)
    seed : int, optional
        Random seed

    Returns
    -------
    x_samples : np.ndarray, shape (n_samples, 1)
        Power samples
    p_samples : np.ndarray, shape (n_samples, n_params)
        Parameter samples
    """
    if lhs is None:
        raise ImportError("pyDOE is required: pip install pyDOE")

    if seed is not None:
        np.random.seed(seed)

    n_params = len(param_bounds)
    total_dim = 1 + n_params  # power + parameters

    # LHS sampling in [0, 1]
    samples = lhs(total_dim, n_samples)

    # Scale power
    power_min, power_max = power_range
    x_samples = samples[:, 0:1] * (power_max - power_min) + power_min

    # Scale parameters
    p_samples = np.zeros((n_samples, n_params))
    for i, (p_min, p_max) in enumerate(param_bounds):
        p_samples[:, i] = samples[:, i + 1] * (p_max - p_min) + p_min

    return x_samples, p_samples


def filter_valid_samples(data: np.ndarray) -> np.ndarray:
    """
    Filter out invalid samples containing NaN

    Parameters
    ----------
    data : np.ndarray
        Raw data

    Returns
    -------
    np.ndarray
        Filtered valid data
    """
    valid_mask = ~np.isnan(data).any(axis=1)
    n_valid = valid_mask.sum()
    n_total = len(data)

    if n_valid < n_total:
        print(f"Filtered invalid samples: {n_total - n_valid}/{n_total} samples removed")

    return data[valid_mask]
