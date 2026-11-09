"""
Objective function module

Defines residual computation, RMSE, and the Bayes unified NRMSE percentage objective.
"""

from typing import Optional, Tuple, Union

import numpy as np

# Penalty value used when simulation fails
PENALTY_VALUE = 1e8


def _prepare_residuals(predictions: np.ndarray, observations: np.ndarray) -> np.ndarray:
    """Compute residuals and replace NaN values with the penalty value."""
    residuals = np.asarray(predictions, dtype=float) - np.asarray(observations, dtype=float)
    nan_mask = np.isnan(residuals)
    if np.any(nan_mask):
        residuals[nan_mask] = PENALTY_VALUE
    return residuals


def _compute_bayes_scales(observations: np.ndarray) -> np.ndarray:
    """Compute the normalization scale for each output according to the Bayes objective rule."""
    obs = np.asarray(observations, dtype=float)
    if obs.ndim == 1:
        obs = obs.reshape(1, -1)
    scales = np.mean(np.abs(obs), axis=0)
    scales = np.where(scales < 1e-10, 1.0, scales)
    return scales


def compute_bayes_normalized_residuals(
    predictions: np.ndarray,
    observations: np.ndarray,
    weights: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Compute the normalized residual matrix corresponding to the Bayes unified objective.

    If weights are provided, only output dimensions with weight > 0 are retained;
    if all weights are 0, all dimensions are retained.
    """
    residuals = _prepare_residuals(predictions, observations)
    scales = _compute_bayes_scales(observations)
    normalized_residuals = residuals / scales

    if weights is not None:
        weights_arr = np.asarray(weights, dtype=float).reshape(-1)
        active_mask = weights_arr > 0
        if np.any(active_mask):
            normalized_residuals = normalized_residuals[:, active_mask]

    return normalized_residuals


def compute_bayes_nrmse_percent(
    predictions: np.ndarray,
    observations: np.ndarray,
    weights: Optional[np.ndarray] = None,
) -> Tuple[float, np.ndarray]:
    """
    Compute the Bayes unified objective: normalized RMSE percentage.

    Returns
    -------
    tuple
        (overall_cost_percent, per_output_percent)
    """
    residuals = _prepare_residuals(predictions, observations)
    scales = _compute_bayes_scales(observations)
    normalized_residuals = residuals / scales

    per_output = np.sqrt(np.mean(normalized_residuals**2, axis=0)) * 100.0

    if weights is not None:
        weights_arr = np.asarray(weights, dtype=float).reshape(-1)
        active_mask = weights_arr > 0
        if np.any(active_mask):
            aggregate = normalized_residuals[:, active_mask]
        else:
            aggregate = normalized_residuals
    else:
        aggregate = normalized_residuals

    cost = float(np.sqrt(np.mean(aggregate**2)) * 100.0)
    return cost, per_output


def compute_residuals(
    predictions: np.ndarray,
    observations: np.ndarray,
    weights: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Compute residuals

    Parameters
    ----------
    predictions : np.ndarray, shape (n, m)
        Model predictions
    observations : np.ndarray, shape (n, m)
        Experimental observations
    weights : np.ndarray, optional
        Weight matrix with the same shape as predictions

    Returns
    -------
    residuals : np.ndarray
        Flattened weighted residual vector
    """
    residuals = _prepare_residuals(predictions, observations)

    if weights is not None:
        residuals = residuals * weights

    return residuals.flatten()


def compute_nrmse_percent(
    predictions: np.ndarray,
    observations: np.ndarray,
    weights: Optional[np.ndarray] = None,
) -> float:
    """
    Compute the normalized root mean square error percentage (NRMSE%)

    Parameters
    ----------
    predictions : np.ndarray, shape (n_points, 3)
        Model predictions [width, depth, area]
    observations : np.ndarray, shape (n_points, 3)
        Experimental observations
    weights : np.ndarray, shape (3,), optional
        Output weights [w_width, w_depth, w_area]
        e.g. [1.0, 1.0, 0.0] means area is not calibrated

    Returns
    -------
    float
        NRMSE percentage
    """
    cost, _ = compute_bayes_nrmse_percent(
        predictions,
        observations,
        weights=weights,
    )
    return cost


def compute_rmse(
    predictions: np.ndarray,
    observations: np.ndarray,
) -> Union[float, np.ndarray]:
    """
    Compute the Root Mean Square Error (RMSE)

    Parameters
    ----------
    predictions : np.ndarray, shape (n, m)
        Model predictions
    observations : np.ndarray, shape (n, m)
        Experimental observations

    Returns
    -------
    float or np.ndarray
        Scalar if input is 1D; otherwise returns per-column RMSE array
    """
    residuals = _prepare_residuals(predictions, observations)

    if residuals.ndim == 1:
        return float(np.sqrt(np.mean(residuals**2)))
    return np.sqrt(np.mean(residuals**2, axis=0))


def compute_mae(
    predictions: np.ndarray,
    observations: np.ndarray,
) -> Union[float, np.ndarray]:
    """
    Compute the Mean Absolute Error (MAE)

    Parameters
    ----------
    predictions : np.ndarray
        Model predictions
    observations : np.ndarray
        Experimental observations

    Returns
    -------
    float or np.ndarray
        MAE value
    """
    residuals = np.abs(np.asarray(predictions, dtype=float) - np.asarray(observations, dtype=float))

    residuals = np.where(np.isnan(residuals), PENALTY_VALUE, residuals)

    if residuals.ndim == 1:
        return float(np.mean(residuals))
    return np.mean(residuals, axis=0)


def compute_relative_error(
    predictions: np.ndarray,
    observations: np.ndarray,
) -> np.ndarray:
    """
    Compute the relative error (%)

    Parameters
    ----------
    predictions : np.ndarray
        Model predictions
    observations : np.ndarray
        Experimental observations

    Returns
    -------
    np.ndarray
        Relative error (%)
    """
    denom = np.abs(observations) + 1e-10
    return 100 * np.abs(predictions - observations) / denom


def normalized_residuals(
    predictions: np.ndarray,
    observations: np.ndarray,
    scales: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Compute normalized residuals

    Parameters
    ----------
    predictions : np.ndarray
        Model predictions
    observations : np.ndarray
        Experimental observations
    scales : np.ndarray, optional
        Scaling factor for each output; defaults to the standard deviation of observations

    Returns
    -------
    np.ndarray
        Normalized residuals
    """
    residuals = np.asarray(predictions, dtype=float) - np.asarray(observations, dtype=float)

    if scales is None:
        scales = np.std(observations, axis=0)
        scales = np.where(scales < 1e-10, 1.0, scales)

    return residuals / scales
