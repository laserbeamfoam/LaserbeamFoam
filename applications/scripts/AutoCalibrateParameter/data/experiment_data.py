"""
实验数据加载模块
"""

from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd


def load_experiment_data(filepath: Path) -> np.ndarray:
    """
    加载实验数据

    Parameters
    ----------
    filepath : Path
        实验数据 CSV 文件路径

    Returns
    -------
    np.ndarray, shape (n_points, 4)
        实验数据 [power_W, width_um, depth_um, area_um2]
        注意：如果CSV中缺少area列，将自动填充为0
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"实验数据文件不存在: {filepath}")

    df = pd.read_csv(filepath)

    # 检查必需列（power, width, depth是必需的，area是可选的）
    required_cols = ["power_W", "width_um", "depth_um"]
    optional_cols = ["area_um2"]

    # 尝试列名映射
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

    # 检查必需列
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        raise ValueError(f"实验数据缺少必需列: {missing_cols}")

    # 如果缺少area列，填充为0
    if "area_um2" not in df.columns:
        print("警告: 实验数据中缺少 area_um2 列，将填充为0")
        df["area_um2"] = 0.0

    return df[["power_W", "width_um", "depth_um", "area_um2"]].values


def prepare_multivariate_data(filepath: Path) -> np.ndarray:
    """
    准备多变量实验数据

    Parameters
    ----------
    filepath : Path
        实验数据 CSV 文件路径

    Returns
    -------
    np.ndarray, shape (n_points, 4)
        实验数据 [power, width, depth, area]
    """
    print(f"读取实验数据（多变量格式）: {filepath}")
    data = load_experiment_data(filepath)
    print(f"实验数据形状: {data.shape}, 功率点数: {len(data)}\n")
    return data


def split_experiment_data(
    data: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    分离实验数据为输入和输出

    Parameters
    ----------
    data : np.ndarray, shape (n, 4)
        完整实验数据

    Returns
    -------
    x : np.ndarray, shape (n, 1)
        输入 (功率)
    y : np.ndarray, shape (n, 3)
        输出 (宽度, 深度, 面积)
    """
    x = data[:, 0:1]  # 功率
    y = data[:, 1:4]  # 宽度、深度、面积
    return x, y


def get_power_points(data: np.ndarray) -> np.ndarray:
    """
    获取实验数据中的功率点列表

    Parameters
    ----------
    data : np.ndarray
        实验数据

    Returns
    -------
    np.ndarray
        唯一功率点列表
    """
    return np.unique(data[:, 0])
