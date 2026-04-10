"""
熔池仿真参数定义

参数元数据统一来自参数注册中心（utils.param_registry）。
"""

from typing import Dict, List, Tuple

from ..utils.param_registry import (
    DEFAULT_BOUNDS,
    PARAM_LABELS,
    PARAM_LATEX_LABELS,
    PARAM_NAMES,
)

# 输出名称
OUTPUT_NAMES: List[str] = ["width", "depth"]

# 输出标签 (用于绘图)
OUTPUT_LABELS: Dict[str, str] = {
    "width": r"Width ($\mu$m)",
    "depth": r"Depth ($\mu$m)",
    "area": r"Area ($\mu$m$^2$)",
}

# 输出单位换算因子 (从 SI 到显示单位)
OUTPUT_SCALE_FACTORS: Dict[str, float] = {
    "width": 1e6,   # m -> μm
    "depth": 1e6,   # m -> μm
    "area": 1e12,   # m² -> μm²
}


def get_param_bounds_list() -> List[Tuple[float, float]]:
    """
    返回参数边界列表（按 PARAM_NAMES 顺序）

    Returns
    -------
    list of tuple
        [(sigma_min, sigma_max), (marangoni_min, marangoni_max), ...]
    """
    return [DEFAULT_BOUNDS[name] for name in PARAM_NAMES]


def param_index(name: str) -> int:
    """
    获取参数索引

    Parameters
    ----------
    name : str
        参数名称

    Returns
    -------
    int
        参数在数组中的索引
    """
    return PARAM_NAMES.index(name)
