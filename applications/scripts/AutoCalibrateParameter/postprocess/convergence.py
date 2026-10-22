"""
Convergence plot drawing
"""

from pathlib import Path
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np

from .plots import setup_plot_style, save_figure


def save_convergence_plot(
    output_dir: Path,
    costs: List[float],
    n_initial: int = 0,
    title: str = "Optimization Convergence",
) -> Path:
    """
    Save convergence curve plot

    Parameters
    ----------
    output_dir : Path
        Output directory
    costs : list
        Objective function value at each iteration
    n_initial : int
        Number of initial sampling points (for split display)
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

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    iterations = np.arange(1, len(costs) + 1)
    costs_array = np.array(costs)

    # Cumulative best
    best_so_far = np.minimum.accumulate(costs_array)

    # Left plot: all evaluations
    ax1.semilogy(iterations, costs_array, "o-", alpha=0.5, markersize=3, label="NRMSE")
    ax1.semilogy(iterations, best_so_far, "r-", lw=2, label="Best so far")

    if n_initial > 0 and n_initial < len(costs):
        ax1.axvline(n_initial, color="gray", linestyle="--", alpha=0.7, label="End of initial sampling")

    ax1.set_xlabel("Evaluation")
    ax1.set_ylabel("NRMSE (log scale)")
    ax1.set_title("All Evaluations")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Right plot: best value evolution
    ax2.semilogy(iterations, best_so_far, "r-", lw=2)
    ax2.fill_between(iterations, best_so_far, alpha=0.3)
    ax2.set_xlabel("Evaluation")
    ax2.set_ylabel("Best NRMSE (log scale)")
    ax2.set_title("Best NRMSE Evolution")
    ax2.grid(True, alpha=0.3)

    # Add statistics
    textstr = f"Final best: {best_so_far[-1]:.4e}\nTotal evals: {len(costs)}"
    ax2.text(0.95, 0.95, textstr, transform=ax2.transAxes, fontsize=10,
             verticalalignment="top", horizontalalignment="right",
             bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

    fig.suptitle(title)
    fig.tight_layout()

    output_path = output_dir / "convergence"
    save_figure(fig, output_path)
    plt.close(fig)

    return output_path.with_suffix(".png")


def plot_parameter_evolution(
    output_dir: Path,
    params_history: List[np.ndarray],
    param_names: List[str],
    title: str = "Parameter Evolution",
) -> Path:
    """
    Plot parameter evolution chart

    Parameters
    ----------
    output_dir : Path
        Output directory
    params_history : list of np.ndarray
        Parameter history records
    param_names : list
        Parameter names
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

    params_array = np.array(params_history)
    n_iters, n_params = params_array.shape

    fig, axes = plt.subplots(n_params, 1, figsize=(10, 3 * n_params), sharex=True)
    if n_params == 1:
        axes = [axes]

    iterations = np.arange(1, n_iters + 1)

    for i, (ax, name) in enumerate(zip(axes, param_names)):
        ax.plot(iterations, params_array[:, i], "b-", alpha=0.7)
        ax.scatter(iterations, params_array[:, i], c=iterations, cmap="viridis", s=10)

        # Mark the final value
        ax.axhline(params_array[-1, i], color="r", linestyle="--", alpha=0.5)

        ax.set_ylabel(name)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Iteration")
    fig.suptitle(title)
    fig.tight_layout()

    output_path = output_dir / "parameter_evolution"
    save_figure(fig, output_path)
    plt.close(fig)

    return output_path.with_suffix(".png")


def plot_multi_direction_comparison(
    output_dir: Path,
    results: List[dict],
    title: str = "Multi-Direction Optimization Comparison",
) -> Path:
    """
    Plot multi-direction optimization comparison chart

    Parameters
    ----------
    output_dir : Path
        Output directory
    results : list of dict
        Optimization results for each direction
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

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Left plot: final NRMSE per direction
    direction_ids = [r["direction_id"] for r in results]
    final_costs = [r["cost"] for r in results]

    colors = ["green" if c == min(final_costs) else "steelblue" for c in final_costs]
    ax1.bar(direction_ids, final_costs, color=colors, edgecolor="black")
    ax1.set_xlabel("Direction ID")
    ax1.set_ylabel("Final NRMSE")
    ax1.set_title("Final NRMSE by Direction")
    ax1.set_yscale("log")

    # Mark the best
    best_idx = np.argmin(final_costs)
    ax1.annotate(
        f"Best: {final_costs[best_idx]:.4e}",
        xy=(direction_ids[best_idx], final_costs[best_idx]),
        xytext=(10, 10),
        textcoords="offset points",
        fontsize=10,
        arrowprops=dict(arrowstyle="->", color="red"),
    )

    # Right plot: convergence curves per direction
    for r in results:
        if "cost_history" in r:
            history = r["cost_history"]
            label = f"Dir {r['direction_id']}"
            if r["direction_id"] == direction_ids[best_idx]:
                ax2.semilogy(history, "r-", lw=2, label=label + " (best)")
            else:
                ax2.semilogy(history, alpha=0.5, lw=1, label=label)

    ax2.set_xlabel("Evaluation")
    ax2.set_ylabel("NRMSE (log scale)")
    ax2.set_title("Convergence by Direction")
    ax2.legend(fontsize=8, ncol=2)
    ax2.grid(True, alpha=0.3)

    fig.suptitle(title)
    fig.tight_layout()

    output_path = output_dir / "multi_direction_comparison"
    save_figure(fig, output_path)
    plt.close(fig)

    return output_path.with_suffix(".png")
