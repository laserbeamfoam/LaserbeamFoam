"""
Gradient optimization configuration

Parameter configuration for scipy.optimize.least_squares point estimation optimization
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from .base_config import BaseConfig


@dataclass
class GradientConfig(BaseConfig):
    """
    Gradient optimization configuration (scipy.optimize.least_squares)

    Inherits BaseConfig and adds gradient optimizer parameters

    Attributes
    ----------
    max_nfev : int
        Maximum number of function evaluations
    ftol : float
        Relative convergence tolerance for the objective function
    xtol : float
        Relative convergence tolerance for parameters
    diff_step : float
        Finite difference step size
    k_lower : float
        Lower bound of parameter scaling factor
    k_upper : float
        Upper bound of parameter scaling factor
    n_directions : int
        Number of directions for multi-start optimization
    init_strategy : str
        Initialization strategy
    direction_scale : float
        Direction scaling factor
    max_workers : int
        Maximum number of parallel optimization processes
    log_every_eval : int
        Log interval in number of evaluations
    monitor_interval_sec : float
        Monitor thread polling interval (seconds)
    name : str
        Optimizer name
    """

    # === Optimization parameters ===
    max_nfev: int = 300
    ftol: float = 1e-4
    xtol: float = 1e-4
    diff_step: float = 0.05

    # === Parameter scaling factor bounds ===
    k_lower: float = 0.5
    k_upper: float = 1.5

    # === Multi-start optimization ===
    n_directions: int = 8
    init_strategy: Literal["legacy_single", "legacy_multi", "directional"] = "legacy_multi"
    direction_scale: float = 0.5

    # === Parallelism ===
    max_workers: int = 8

    # === Monitoring ===
    log_every_eval: int = 1
    monitor_interval_sec: float = 0.5

    # === Output ===
    name: str = "meltpool_gradient"

    def to_dict(self) -> dict:
        """Convert to dictionary"""
        base_dict = super().to_dict()
        base_dict.update({
            "max_nfev": self.max_nfev,
            "ftol": self.ftol,
            "xtol": self.xtol,
            "diff_step": self.diff_step,
            "k_lower": self.k_lower,
            "k_upper": self.k_upper,
            "n_directions": self.n_directions,
            "init_strategy": self.init_strategy,
            "direction_scale": self.direction_scale,
            "max_workers": self.max_workers,
            "log_every_eval": self.log_every_eval,
            "monitor_interval_sec": self.monitor_interval_sec,
            "name": self.name,
        })
        return base_dict
