"""
Parameter filtering utilities

Supports selective parameter optimization
"""
from typing import Dict, List, Optional, Tuple

import numpy as np

from .param_registry import canonical_param_name


def filter_params_for_optimization(
    all_param_names: List[str],
    all_bounds: Dict[str, Tuple[float, float]],
    active_params: Optional[List[str]] = None,
    fixed_values: Optional[Dict[str, float]] = None,
) -> Tuple[List[str], List[Tuple[float, float]], Dict[str, float]]:
    """
    Filter parameters to support selective optimization

    Parameters
    ----------
    all_param_names : list of str
        All parameter names (in order)
    all_bounds : dict
        Bounds for all parameters
    active_params : list of str, optional
        Parameter names to optimize; None means optimize all parameters
    fixed_values : dict, optional
        Values for fixed parameters

    Returns
    -------
    active_names : list of str
        Parameter names to optimize
    active_bounds : list of tuple
        Bounds for parameters to optimize
    fixed_dict : dict
        Mapping of fixed parameter names to values
    """
    if active_params is None:
        # Default: optimize all parameters
        active_names = all_param_names.copy()
    else:
        active_names = []
        for name in active_params:
            canonical_name = canonical_param_name(name)
            if canonical_name not in active_names:
                active_names.append(canonical_name)

    # Build optimization space
    active_bounds = [all_bounds[name] for name in active_names]

    # Build fixed parameter dictionary
    fixed_dict = {}
    if fixed_values:
        for name, value in fixed_values.items():
            fixed_dict[canonical_param_name(name)] = value

    # For parameters not in active list, use midpoint of bounds if no fixed value is specified
    for name in all_param_names:
        if name not in active_names and name not in fixed_dict:
            bounds = all_bounds[name]
            fixed_dict[name] = (bounds[0] + bounds[1]) / 2.0

    return active_names, active_bounds, fixed_dict


def merge_active_and_fixed_params(
    all_param_names: List[str],
    active_names: List[str],
    active_values: np.ndarray,
    fixed_dict: Dict[str, float],
) -> np.ndarray:
    """
    Merge active parameter values and fixed parameter values into a complete parameter vector

    Parameters
    ----------
    all_param_names : list of str
        All parameter names (in order)
    active_names : list of str
        Parameter names to optimize
    active_values : np.ndarray
        Values of optimized parameters (corresponding to active_names)
    fixed_dict : dict
        Mapping of fixed parameter names to values

    Returns
    -------
    np.ndarray
        Complete parameter vector (in all_param_names order)
    """
    full_params = {}

    # Fill active parameters
    for name, val in zip(active_names, active_values):
        full_params[name] = val

    # Fill fixed parameters
    for name, val in fixed_dict.items():
        full_params[name] = val

    # Return in all_param_names order
    return np.array([full_params[name] for name in all_param_names])
