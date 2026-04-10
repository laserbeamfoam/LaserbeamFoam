"""
参数过滤工具

用于支持选择性参数优化
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
    过滤参数以支持选择性优化

    Parameters
    ----------
    all_param_names : list of str
        所有参数名称（按顺序）
    all_bounds : dict
        所有参数的边界
    active_params : list of str, optional
        要优化的参数名称，None表示优化所有参数
    fixed_values : dict, optional
        固定参数的值

    Returns
    -------
    active_names : list of str
        要优化的参数名称
    active_bounds : list of tuple
        要优化的参数边界
    fixed_dict : dict
        固定参数的名称->值映射
    """
    if active_params is None:
        # 默认优化所有参数
        active_names = all_param_names.copy()
    else:
        active_names = []
        for name in active_params:
            canonical_name = canonical_param_name(name)
            if canonical_name not in active_names:
                active_names.append(canonical_name)

    # 构建优化空间
    active_bounds = [all_bounds[name] for name in active_names]

    # 构建固定参数字典
    fixed_dict = {}
    if fixed_values:
        for name, value in fixed_values.items():
            fixed_dict[canonical_param_name(name)] = value

    # 对于未在active中的参数，如果没有指定固定值，使用边界中点
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
    合并active参数值和fixed参数值为完整参数向量

    Parameters
    ----------
    all_param_names : list of str
        所有参数名称（按顺序）
    active_names : list of str
        要优化的参数名称
    active_values : np.ndarray
        优化参数的值（与active_names对应）
    fixed_dict : dict
        固定参数的名称->值映射

    Returns
    -------
    np.ndarray
        完整参数向量（按all_param_names顺序）
    """
    full_params = {}

    # 填充active参数
    for name, val in zip(active_names, active_values):
        full_params[name] = val

    # 填充fixed参数
    for name, val in fixed_dict.items():
        full_params[name] = val

    # 按all_param_names顺序返回
    return np.array([full_params[name] for name in all_param_names])
