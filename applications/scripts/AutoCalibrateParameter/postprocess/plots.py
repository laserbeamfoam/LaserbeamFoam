"""
通用绑图工具
"""

from pathlib import Path
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np


def setup_plot_style() -> None:
    """设置绘图风格"""
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update({
        "figure.figsize": (10, 6),
        "figure.dpi": 100,
        "font.size": 12,
        "axes.labelsize": 12,
        "axes.titlesize": 14,
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
    })


def save_figure(
    fig: plt.Figure,
    filepath: Path,
    formats: List[str] = ["png"],
    dpi: int = 150,
) -> None:
    """
    保存图形到多种格式

    Parameters
    ----------
    fig : matplotlib.Figure
        图形对象
    filepath : Path
        输出路径（不含扩展名）
    formats : list
        输出格式列表
    dpi : int
        分辨率
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    for fmt in formats:
        output_path = filepath.with_suffix(f".{fmt}")
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        print(f"  已保存: {output_path}")


def plot_parameter_trace(
    samples: np.ndarray,
    param_names: List[str],
    output_path: Optional[Path] = None,
    title: str = "Parameter Trace",
) -> plt.Figure:
    """
    绘制参数 trace 图

    Parameters
    ----------
    samples : np.ndarray, shape (n_samples, n_params)
        参数样本
    param_names : list
        参数名称列表
    output_path : Path, optional
        输出路径
    title : str
        标题

    Returns
    -------
    matplotlib.Figure
    """
    setup_plot_style()

    n_params = samples.shape[1]
    fig, axes = plt.subplots(n_params, 1, figsize=(10, 2.5 * n_params), sharex=True)

    if n_params == 1:
        axes = [axes]

    for i, (ax, name) in enumerate(zip(axes, param_names)):
        ax.plot(samples[:, i], alpha=0.7, lw=0.5)
        ax.set_ylabel(name)
        ax.axhline(np.mean(samples[:, i]), color="r", linestyle="--", lw=1, label="Mean")

    axes[-1].set_xlabel("Iteration")
    fig.suptitle(title)
    fig.tight_layout()

    if output_path:
        save_figure(fig, output_path)

    return fig


def plot_parameter_histogram(
    samples: np.ndarray,
    param_names: List[str],
    output_path: Optional[Path] = None,
    title: str = "Parameter Distribution",
    bins: int = 50,
) -> plt.Figure:
    """
    绘制参数直方图

    Parameters
    ----------
    samples : np.ndarray, shape (n_samples, n_params)
        参数样本
    param_names : list
        参数名称列表
    output_path : Path, optional
        输出路径
    title : str
        标题
    bins : int
        直方图 bin 数量

    Returns
    -------
    matplotlib.Figure
    """
    setup_plot_style()

    n_params = samples.shape[1]
    n_cols = min(2, n_params)
    n_rows = (n_params + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    axes = np.atleast_1d(axes).flatten()

    for i, (ax, name) in enumerate(zip(axes[:n_params], param_names)):
        data = samples[:, i]
        ax.hist(data, bins=bins, density=True, alpha=0.7, edgecolor="black")

        # 统计信息
        mean = np.mean(data)
        median = np.median(data)
        std = np.std(data)
        q05, q95 = np.percentile(data, [5, 95])

        ax.axvline(mean, color="r", linestyle="-", lw=2, label=f"Mean: {mean:.4g}")
        ax.axvline(median, color="g", linestyle="--", lw=2, label=f"Median: {median:.4g}")
        ax.axvspan(q05, q95, alpha=0.2, color="blue", label=f"90% CI")

        ax.set_xlabel(name)
        ax.set_ylabel("Density")
        ax.legend(fontsize=8)

    # 隐藏多余的子图
    for ax in axes[n_params:]:
        ax.set_visible(False)

    fig.suptitle(title)
    fig.tight_layout()

    if output_path:
        save_figure(fig, output_path)

    return fig
