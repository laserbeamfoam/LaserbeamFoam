"""
合成数据管理模块
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
    加载合成数据

    Parameters
    ----------
    filepath : Path
        合成数据文件路径

    Returns
    -------
    np.ndarray, shape (n_samples, 9)
        [power, sigma, marangoni, substrate_temp, absorptivity, recoilCoeff, radius_flavour, width, depth, area]
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"合成数据文件不存在: {filepath}")

    data = np.loadtxt(filepath)
    print(f"加载合成数据: {filepath}, 形状 {data.shape}")
    return data


def save_synthetic_data(
    filepath: Path,
    data: np.ndarray,
    header: str = "power sigma marangoni substrate_temp absorptivity recoilCoeff radius_flavour width depth area",
) -> None:
    """
    保存合成数据

    Parameters
    ----------
    filepath : Path
        输出文件路径
    data : np.ndarray
        合成数据
    header : str
        文件头注释
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(filepath, data, header=header)
    print(f"合成数据已保存至: {filepath}, 形状: {data.shape}")


def append_synthetic_data(filepath: Path, new_data: np.ndarray) -> np.ndarray:
    """
    追加新数据到合成数据文件

    Parameters
    ----------
    filepath : Path
        合成数据文件路径
    new_data : np.ndarray
        新数据

    Returns
    -------
    np.ndarray
        更新后的完整数据
    """
    import shutil

    existing_data = load_synthetic_data(filepath)
    updated_data = np.vstack([existing_data, new_data])

    # 备份原文件
    backup_path = filepath.with_suffix(filepath.suffix + ".backup")
    shutil.copy(filepath, backup_path)

    # 保存更新后的数据
    save_synthetic_data(filepath, updated_data)
    print(f"数据已更新: {filepath}, 新形状 {updated_data.shape}")

    return updated_data


def generate_lhs_samples(
    n_samples: int,
    param_bounds: List[Tuple[float, float]],
    power_range: Tuple[float, float] = (140.0, 260.0),
    seed: Optional[int] = 42,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    使用拉丁超立方采样生成样本

    Parameters
    ----------
    n_samples : int
        样本数量
    param_bounds : list of tuple
        参数边界列表 [(min, max), ...]
    power_range : tuple
        功率范围 (min, max)
    seed : int, optional
        随机种子

    Returns
    -------
    x_samples : np.ndarray, shape (n_samples, 1)
        功率样本
    p_samples : np.ndarray, shape (n_samples, n_params)
        参数样本
    """
    if lhs is None:
        raise ImportError("需要安装 pyDOE: pip install pyDOE")

    if seed is not None:
        np.random.seed(seed)

    n_params = len(param_bounds)
    total_dim = 1 + n_params  # 功率 + 参数

    # LHS 采样 [0, 1]
    samples = lhs(total_dim, n_samples)

    # 缩放功率
    power_min, power_max = power_range
    x_samples = samples[:, 0:1] * (power_max - power_min) + power_min

    # 缩放参数
    p_samples = np.zeros((n_samples, n_params))
    for i, (p_min, p_max) in enumerate(param_bounds):
        p_samples[:, i] = samples[:, i + 1] * (p_max - p_min) + p_min

    return x_samples, p_samples


def filter_valid_samples(data: np.ndarray) -> np.ndarray:
    """
    过滤掉包含 NaN 的无效样本

    Parameters
    ----------
    data : np.ndarray
        原始数据

    Returns
    -------
    np.ndarray
        过滤后的有效数据
    """
    valid_mask = ~np.isnan(data).any(axis=1)
    n_valid = valid_mask.sum()
    n_total = len(data)

    if n_valid < n_total:
        print(f"过滤无效样本: {n_total - n_valid}/{n_total} 个样本被移除")

    return data[valid_mask]
