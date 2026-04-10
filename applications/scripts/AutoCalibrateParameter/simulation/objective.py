"""
目标函数模块

定义残差计算、RMSE 及 Bayes 统一 NRMSE 百分比目标。
"""

from typing import Optional, Tuple, Union

import numpy as np

# 仿真失败时的惩罚值
PENALTY_VALUE = 1e8


def _prepare_residuals(predictions: np.ndarray, observations: np.ndarray) -> np.ndarray:
    """计算残差并将 NaN 替换为惩罚值。"""
    residuals = np.asarray(predictions, dtype=float) - np.asarray(observations, dtype=float)
    nan_mask = np.isnan(residuals)
    if np.any(nan_mask):
        residuals[nan_mask] = PENALTY_VALUE
    return residuals


def _compute_bayes_scales(observations: np.ndarray) -> np.ndarray:
    """按 Bayes 目标规则计算每个输出的归一化尺度。"""
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
    计算 Bayes 统一目标对应的归一化残差矩阵。

    若给定 weights，则仅保留权重大于 0 的输出维度；若全部为 0，则保留全部维度。
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
    计算 Bayes 统一目标：归一化 RMSE 百分比。

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
    计算残差

    Parameters
    ----------
    predictions : np.ndarray, shape (n, m)
        模型预测值
    observations : np.ndarray, shape (n, m)
        实验观测值
    weights : np.ndarray, optional
        权重矩阵，与 predictions 形状相同

    Returns
    -------
    residuals : np.ndarray
        展平后的加权残差向量
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
    计算归一化均方根误差百分比 (NRMSE%)

    Parameters
    ----------
    predictions : np.ndarray, shape (n_points, 3)
        模型预测值 [width, depth, area]
    observations : np.ndarray, shape (n_points, 3)
        实验观测值
    weights : np.ndarray, shape (3,), optional
        输出权重 [w_width, w_depth, w_area]
        例如 [1.0, 1.0, 0.0] 表示不校准area

    Returns
    -------
    float
        NRMSE 百分比
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
    计算均方根误差 (Root Mean Square Error)

    Parameters
    ----------
    predictions : np.ndarray, shape (n, m)
        模型预测值
    observations : np.ndarray, shape (n, m)
        实验观测值

    Returns
    -------
    float or np.ndarray
        如果输入是 1D，返回标量；否则按列返回 RMSE 数组
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
    计算平均绝对误差 (Mean Absolute Error)

    Parameters
    ----------
    predictions : np.ndarray
        模型预测值
    observations : np.ndarray
        实验观测值

    Returns
    -------
    float or np.ndarray
        MAE 值
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
    计算相对误差 (%)

    Parameters
    ----------
    predictions : np.ndarray
        模型预测值
    observations : np.ndarray
        实验观测值

    Returns
    -------
    np.ndarray
        相对误差 (%)
    """
    denom = np.abs(observations) + 1e-10
    return 100 * np.abs(predictions - observations) / denom


def normalized_residuals(
    predictions: np.ndarray,
    observations: np.ndarray,
    scales: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    计算归一化残差

    Parameters
    ----------
    predictions : np.ndarray
        模型预测值
    observations : np.ndarray
        实验观测值
    scales : np.ndarray, optional
        各输出的缩放因子，默认使用观测值的标准差

    Returns
    -------
    np.ndarray
        归一化残差
    """
    residuals = np.asarray(predictions, dtype=float) - np.asarray(observations, dtype=float)

    if scales is None:
        scales = np.std(observations, axis=0)
        scales = np.where(scales < 1e-10, 1.0, scales)

    return residuals / scales
