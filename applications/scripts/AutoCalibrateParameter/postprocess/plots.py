"""
General plotting utilities
"""

from pathlib import Path
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np


def setup_plot_style() -> None:
    """Set up the plot style"""
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
    Save figure in multiple formats

    Parameters
    ----------
    fig : matplotlib.Figure
        Figure object
    filepath : Path
        Output path (without extension)
    formats : list
        List of output formats
    dpi : int
        Resolution
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    for fmt in formats:
        output_path = filepath.with_suffix(f".{fmt}")
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
        print(f"  Saved: {output_path}")


def plot_parameter_trace(
    samples: np.ndarray,
    param_names: List[str],
    output_path: Optional[Path] = None,
    title: str = "Parameter Trace",
) -> plt.Figure:
    """
    Plot parameter trace chart

    Parameters
    ----------
    samples : np.ndarray, shape (n_samples, n_params)
        Parameter samples
    param_names : list
        List of parameter names
    output_path : Path, optional
        Output path
    title : str
        Title

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
    Plot parameter histogram

    Parameters
    ----------
    samples : np.ndarray, shape (n_samples, n_params)
        Parameter samples
    param_names : list
        List of parameter names
    output_path : Path, optional
        Output path
    title : str
        Title
    bins : int
        Number of histogram bins

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

        # Statistics
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

    # Hide extra subplots
    for ax in axes[n_params:]:
        ax.set_visible(False)

    fig.suptitle(title)
    fig.tight_layout()

    if output_path:
        save_figure(fig, output_path)

    return fig
