"""
Bayesian optimization configuration

Parameter configuration for scikit-optimize point estimation optimization
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from .base_config import BaseConfig


@dataclass
class BayesConfig(BaseConfig):
    """
    Bayesian optimization configuration (scikit-optimize)

    Inherits BaseConfig and adds skopt optimizer parameters

    Attributes
    ----------
    n_initial_points : int
        Number of initial random sampling points
    n_batches : int
        Number of Bayesian optimization batches
    batch_size : int
        Number of candidate points evaluated in parallel per batch
    acq_func : str
        Acquisition function ("EI", "PI", "LCB", "gp_hedge")
    xi : float
        Exploration parameter for EI/PI
    kappa : float
        Exploration parameter for LCB
    noise : float
        GP noise parameter
    k_lower : float
        Lower bound of parameter scaling factor
    k_upper : float
        Upper bound of parameter scaling factor
    max_workers : int
        Maximum number of parallel evaluation processes
    resume : bool
        Whether to resume from the last interruption
    name : str
        Optimizer name
    """

    # === Sampling configuration ===
    n_initial_points: int = 15
    n_batches: int = 1
    batch_size: int = 1

    # === Acquisition function ===
    acq_func: Literal["EI", "PI", "LCB", "gp_hedge"] = "EI"
    xi: float = 0.01
    kappa: float = 1.96
    noise: float = 1e-10

    # === Parameter scaling factor bounds ===
    k_lower: float = 0.5
    k_upper: float = 1.5

    # === Parallelism and resume ===
    max_workers: int = 1
    resume: bool = True
    warmstart_csv: str = "bayes_history.csv"  # history data CSV file path for warm start

    # === Output ===
    name: str = "meltpool_bayes"

    def to_dict(self) -> dict:
        """Convert to dictionary"""
        base_dict = super().to_dict()
        base_dict.update({
            "n_initial_points": self.n_initial_points,
            "n_batches": self.n_batches,
            "batch_size": self.batch_size,
            "acq_func": self.acq_func,
            "xi": self.xi,
            "kappa": self.kappa,
            "noise": self.noise,
            "k_lower": self.k_lower,
            "k_upper": self.k_upper,
            "max_workers": self.max_workers,
            "resume": self.resume,
            "warmstart_csv": self.warmstart_csv,
            "name": self.name,
        })
        return base_dict
