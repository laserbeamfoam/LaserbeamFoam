"""
Experiment vs simulation comparison plots
"""

from pathlib import Path
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np

from .plots import setup_plot_style, save_figure


def save_comparison_plot(
    output_dir: Path,
    power_points: np.ndarray,
    exp_observations: np.ndarray,
    sim_predictions: np.ndarray,
    output_names: List[str] = ["Width", "Depth"],
    output_units: List[str] = ["μm", "μm"],
    nrmse: Optional[float] = None,
    title: str = "Experiment vs Simulation",
) -> Path:
    """
    Save experiment vs simulation comparison plot

    Parameters
    ----------
    output_dir : Path
        Output directory
    power_points : np.ndarray
        Power points
    exp_observations : np.ndarray, shape (n, 3)
        Experimental observations [width, depth, area]
    sim_predictions : np.ndarray, shape (n, 3)
        Simulation predictions [width, depth, area]
    output_names : list
        Output names
    output_units : list
        Output units
    nrmse : float, optional
        NRMSE value
    title : str
        Plot title

    Returns
    -------
    Path
        Image path
    """
    setup_plot_style()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_outputs = exp_observations.shape[1]
    fig, axes = plt.subplots(1, n_outputs, figsize=(5 * n_outputs, 4))
    if n_outputs == 1:
        axes = [axes]

    for i, (ax, name, unit) in enumerate(zip(axes, output_names, output_units)):
        exp_vals = exp_observations[:, i]
        sim_vals = sim_predictions[:, i]

        # Plot comparison
        ax.plot(power_points, exp_vals, "ko-", markersize=8, label="Experiment")
        ax.plot(power_points, sim_vals, "rs--", markersize=8, label="Simulation")

        # Fill error region
        ax.fill_between(
            power_points,
            exp_vals,
            sim_vals,
            alpha=0.2,
            color="gray",
        )

        ax.set_xlabel("Power (W)")
        ax.set_ylabel(f"{name} ({unit})")
        ax.legend()
        ax.grid(True, alpha=0.3)

        # Compute per-output RMSE
        rmse = np.sqrt(np.mean((exp_vals - sim_vals) ** 2))
        ax.set_title(f"{name} (RMSE: {rmse:.2f} {unit})")

    if nrmse is not None:
        fig.suptitle(f"{title}\nTotal NRMSE: {nrmse:.4e}")
    else:
        fig.suptitle(title)

    fig.tight_layout()

    output_path = output_dir / "comparison"
    save_figure(fig, output_path)
    plt.close(fig)

    return output_path.with_suffix(".png")


def plot_predictions_vs_experiments(
    output_dir: Path,
    exp_observations: np.ndarray,
    sim_predictions: np.ndarray,
    output_names: List[str] = ["Width", "Depth"],
    output_units: List[str] = ["μm", "μm"],
    title: str = "Predictions vs Experiments",
) -> Path:
    """
    Plot predictions vs experiments scatter plot (45-degree line)

    Parameters
    ----------
    output_dir : Path
        Output directory
    exp_observations : np.ndarray
        Experimental observations
    sim_predictions : np.ndarray
        Simulation predictions
    output_names : list
        Output names
    output_units : list
        Output units
    title : str
        Plot title

    Returns
    -------
    Path
        Image path
    """
    setup_plot_style()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_outputs = exp_observations.shape[1]
    fig, axes = plt.subplots(1, n_outputs, figsize=(5 * n_outputs, 5))
    if n_outputs == 1:
        axes = [axes]

    for i, (ax, name, unit) in enumerate(zip(axes, output_names, output_units)):
        exp_vals = exp_observations[:, i]
        sim_vals = sim_predictions[:, i]

        # Scatter plot
        ax.scatter(exp_vals, sim_vals, s=80, c="steelblue", edgecolor="black", alpha=0.7)

        # 45-degree line
        min_val = min(exp_vals.min(), sim_vals.min())
        max_val = max(exp_vals.max(), sim_vals.max())
        margin = (max_val - min_val) * 0.1
        line_range = [min_val - margin, max_val + margin]
        ax.plot(line_range, line_range, "k--", lw=1, label="Perfect fit")

        # ±10% region
        ax.fill_between(
            line_range,
            [v * 0.9 for v in line_range],
            [v * 1.1 for v in line_range],
            alpha=0.1,
            color="green",
            label="±10%",
        )

        ax.set_xlim(line_range)
        ax.set_ylim(line_range)
        ax.set_xlabel(f"Experiment ({unit})")
        ax.set_ylabel(f"Simulation ({unit})")
        ax.legend()
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.3)

        # R² value
        ss_res = np.sum((exp_vals - sim_vals) ** 2)
        ss_tot = np.sum((exp_vals - np.mean(exp_vals)) ** 2)
        r2 = 1 - ss_res / (ss_tot + 1e-10)
        ax.set_title(f"{name} (R² = {r2:.4f})")

    fig.suptitle(title)
    fig.tight_layout()

    output_path = output_dir / "predictions_vs_experiments"
    save_figure(fig, output_path)
    plt.close(fig)

    return output_path.with_suffix(".png")


def plot_residuals(
    output_dir: Path,
    power_points: np.ndarray,
    residuals: np.ndarray,
    output_names: List[str] = ["Width", "Depth"],
    output_units: List[str] = ["μm", "μm"],
    title: str = "Residuals Analysis",
) -> Path:
    """
    Plot residual analysis chart

    Parameters
    ----------
    output_dir : Path
        Output directory
    power_points : np.ndarray
        Power points
    residuals : np.ndarray
        Residuals (simulation - experiment)
    output_names : list
        Output names
    output_units : list
        Output units
    title : str
        Plot title

    Returns
    -------
    Path
        Image path
    """
    setup_plot_style()

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_outputs = residuals.shape[1]
    fig, axes = plt.subplots(2, n_outputs, figsize=(5 * n_outputs, 8))

    for i, (name, unit) in enumerate(zip(output_names, output_units)):
        res = residuals[:, i]

        # Top plot: residuals vs power
        ax1 = axes[0, i] if n_outputs > 1 else axes[0]
        ax1.bar(power_points, res, width=15, color="steelblue", edgecolor="black")
        ax1.axhline(0, color="r", linestyle="--", lw=1)
        ax1.set_xlabel("Power (W)")
        ax1.set_ylabel(f"Residual ({unit})")
        ax1.set_title(f"{name} Residuals")
        ax1.grid(True, alpha=0.3)

        # Bottom plot: residuals histogram
        ax2 = axes[1, i] if n_outputs > 1 else axes[1]
        ax2.hist(res, bins=10, color="steelblue", edgecolor="black", density=True)
        ax2.axvline(0, color="r", linestyle="--", lw=1)
        ax2.axvline(np.mean(res), color="g", linestyle="-", lw=2, label=f"Mean: {np.mean(res):.2f}")
        ax2.set_xlabel(f"Residual ({unit})")
        ax2.set_ylabel("Density")
        ax2.legend()

    fig.suptitle(title)
    fig.tight_layout()

    output_path = output_dir / "residuals"
    save_figure(fig, output_path)
    plt.close(fig)

    return output_path.with_suffix(".png")
