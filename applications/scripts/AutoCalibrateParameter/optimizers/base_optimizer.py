"""
Optimizer base class

Defines the unified optimizer interface and result data structure
"""

from __future__ import annotations

import json
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np


@dataclass
class OptimizationResult:
    """
    Unified optimization result data structure

    Attributes
    ----------
    method : str
        Optimization method name ("Bayesian", "Gradient")
    best_params : np.ndarray
        Optimal parameters [sigma, marangoni, substrate_temp, absorptivity, recoilCoeff, radius_flavour]
    best_cost : float
        Best objective function value (NRMSE)
    n_evaluations : int
        Total number of evaluations
    history : dict
        Optimization history record
    message : str
        Status message
    extra : dict, optional
        Method-specific extra info
    """

    method: str
    best_params: np.ndarray
    best_cost: float
    n_evaluations: int
    history: Dict[str, Any]
    message: str
    extra: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to a JSON-serializable dictionary"""
        result = {
            "method": self.method,
            "best_params": self.best_params.tolist(),
            "best_cost": self.best_cost,
            "n_evaluations": self.n_evaluations,
            "message": self.message,
        }

        if self.extra:
            result["extra"] = self.extra

        return result

    def save(self, filepath: Path) -> None:
        """Save results to a JSON file"""
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)

        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, indent=2, ensure_ascii=False)

        print(f"Optimization results saved to: {filepath}")

    @classmethod
    def load(cls, filepath: Path) -> "OptimizationResult":
        """Load results from a JSON file"""
        with open(filepath, "r", encoding="utf-8") as f:
            data = json.load(f)

        return cls(
            method=data["method"],
            best_params=np.array(data["best_params"]),
            best_cost=data["best_cost"],
            n_evaluations=data["n_evaluations"],
            history={},
            message=data["message"],
            extra=data.get("extra", {}),
        )


class BaseOptimizer(ABC):
    """
    Abstract optimizer base class

    All optimizers must implement the optimize() method
    """

    def __init__(self, config):
        """
        Initialize the optimizer

        Parameters
        ----------
        config : BaseConfig
            Configuration object
        """
        self.config = config
        self.history: Dict[str, List[Any]] = {
            "params": [],
            "costs": [],
            "iterations": [],
        }

    @abstractmethod
    def get_name(self) -> str:
        """
        Return the optimizer name

        Returns
        -------
        str
            Optimizer name
        """
        pass

    @abstractmethod
    def optimize(self, exp_data: np.ndarray) -> OptimizationResult:
        """
        Execute optimization

        Parameters
        ----------
        exp_data : np.ndarray, shape (n, 4)
            Experimental data [power, width, depth, area]

        Returns
        -------
        OptimizationResult
            Optimization result
        """
        pass

    def save_results(
        self,
        result: OptimizationResult,
        output_dir: Path,
    ) -> None:
        """
        Save optimization results

        Parameters
        ----------
        result : OptimizationResult
            Optimization result
        output_dir : Path
            Output directory
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Save main result
        result.save(output_dir / "opt_result.json")

        # Save history
        if result.history:
            history_file = output_dir / "history.json"
            serializable_history = {}
            for k, v in result.history.items():
                if isinstance(v, np.ndarray):
                    serializable_history[k] = v.tolist()
                elif isinstance(v, list) and len(v) > 0:
                    if isinstance(v[0], np.ndarray):
                        serializable_history[k] = [x.tolist() for x in v]
                    else:
                        serializable_history[k] = v
                else:
                    serializable_history[k] = v

            with open(history_file, "w", encoding="utf-8") as f:
                json.dump(serializable_history, f, indent=2)

    def _log_iteration(
        self,
        params: np.ndarray,
        cost: float,
        iteration: int,
    ) -> None:
        """Record iteration information"""
        self.history["params"].append(params.copy())
        self.history["costs"].append(cost)
        self.history["iterations"].append(iteration)

    def _print_progress(
        self,
        iteration: int,
        total: int,
        cost: float,
        best_cost: float,
    ) -> None:
        """Print progress information"""
        print(
            f"  Iteration {iteration}/{total}: "
            f"current NRMSE = {cost:.4e}, best NRMSE = {best_cost:.4e}"
        )
