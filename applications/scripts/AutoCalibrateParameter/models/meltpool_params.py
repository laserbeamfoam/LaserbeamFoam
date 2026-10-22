"""
Melt pool simulation parameter definitions

Parameter metadata is unified from the parameter registry (utils.param_registry).
"""

from typing import Dict, List, Tuple

from ..utils.param_registry import (
    DEFAULT_BOUNDS,
    PARAM_LABELS,
    PARAM_LATEX_LABELS,
    PARAM_NAMES,
)

# Output names
OUTPUT_NAMES: List[str] = ["width", "depth"]

# Output labels (used for plotting)
OUTPUT_LABELS: Dict[str, str] = {
    "width": r"Width ($\mu$m)",
    "depth": r"Depth ($\mu$m)",
    "area": r"Area ($\mu$m$^2$)",
}

# Output unit conversion factors (from SI to display units)
OUTPUT_SCALE_FACTORS: Dict[str, float] = {
    "width": 1e6,   # m -> μm
    "depth": 1e6,   # m -> μm
    "area": 1e12,   # m² -> μm²
}


def get_param_bounds_list() -> List[Tuple[float, float]]:
    """
    Return the list of parameter bounds (in PARAM_NAMES order)

    Returns
    -------
    list of tuple
        [(sigma_min, sigma_max), (marangoni_min, marangoni_max), ...]
    """
    return [DEFAULT_BOUNDS[name] for name in PARAM_NAMES]


def param_index(name: str) -> int:
    """
    Get the parameter index

    Parameters
    ----------
    name : str
        Parameter name

    Returns
    -------
    int
        Index of the parameter in the array
    """
    return PARAM_NAMES.index(name)
