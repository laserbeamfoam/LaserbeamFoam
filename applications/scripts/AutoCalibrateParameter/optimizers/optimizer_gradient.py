"""
Gradient optimizer

Uses scipy.optimize.least_squares for gradient-based optimization (point estimate)
"""

from __future__ import annotations

import json
import threading
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional

import numpy as np
from scipy.optimize import least_squares

from .base_optimizer import BaseOptimizer, OptimizationResult
from ..simulation.objective import compute_bayes_nrmse_percent, compute_bayes_normalized_residuals

if TYPE_CHECKING:
    from ..config.gradient_config import GradientConfig
    from ..simulation.simulation_runner import SimulationRunner


class GradientOptimizer(BaseOptimizer):
    """
    Gradient optimizer (scipy.optimize.least_squares)

    Features:
    - Gradient estimation based on finite differences
    - Multi-start parallel optimization
    - Returns point estimate

    Attributes
    ----------
    config : GradientConfig
        Gradient optimization configuration
    runner : SimulationRunner
        Simulation runner
    """

    def __init__(self, config: "GradientConfig", runner: "SimulationRunner"):
        """
        Initialize the gradient optimizer

        Parameters
        ----------
        config : GradientConfig
            Gradient optimization configuration
        runner : SimulationRunner
            Simulation runner
        """
        super().__init__(config)
        self.config = config
        self.runner = runner

    def get_name(self) -> str:
        return "Gradient"

    def optimize(self, exp_data: np.ndarray) -> OptimizationResult:
        """
        Execute gradient optimization

        Parameters
        ----------
        exp_data : np.ndarray, shape (n, 4)
            Experimental data [power, width, depth, area]

        Returns
        -------
        OptimizationResult
        """
        cfg = self.config

        print(f"\n{'='*60}")
        print(f"Gradient optimization (scipy.optimize.least_squares)")
        print(f"{'='*60}")
        print(f"Number of directions: {cfg.n_directions}, max evaluations: {cfg.max_nfev}")
        print(f"Initialization strategy: {cfg.init_strategy}")
        print(f"Parallel workers: {cfg.max_workers}")
        print()

        # Run directory
        cfg.runs_root.mkdir(parents=True, exist_ok=True)

        # Store results for all directions
        all_results: List[Dict] = []

        # Start monitoring thread
        stop_event = threading.Event()
        monitor_thread = threading.Thread(
            target=self._monitor_loop,
            args=(stop_event,),
            daemon=True,
        )
        monitor_thread.start()

        try:
            # Run each direction sequentially (simplified, no process pool)
            for direction_id in range(1, cfg.n_directions + 1):
                print(f"\nDirection {direction_id}/{cfg.n_directions}")
                result = self._optimize_single_direction(direction_id, exp_data)
                all_results.append(result)
                print(f"  Done: NRMSE = {result['cost']:.4f}%")

        finally:
            stop_event.set()
            monitor_thread.join(timeout=2.0)

        # Select the best result
        if not all_results:
            return OptimizationResult(
                method="Gradient",
                best_params=np.zeros(len(self.runner.active_names)),
                best_cost=float("inf"),
                n_evaluations=0,
                history={},
                message="All directions failed optimization",
            )

        best = min(all_results, key=lambda r: r["cost"])
        total_evals = sum(r.get("nfev", 0) for r in all_results)

        print(f"\nBest direction: {best['direction_id']}, NRMSE = {best['cost']:.4f}%")

        return OptimizationResult(
            method="Gradient",
            best_params=np.array(best["params"]),
            best_cost=best["cost"],
            n_evaluations=total_evals,
            history={"all_results": all_results},
            message=f"Gradient optimization complete, best direction {best['direction_id']}",
            extra={"all_directions": all_results},
        )

    def _optimize_single_direction(
        self,
        direction_id: int,
        exp_data: np.ndarray,
    ) -> Dict:
        """
        Single-direction optimization

        Parameters
        ----------
        direction_id : int
            Direction ID
        exp_data : np.ndarray
            Experimental data

        Returns
        -------
        dict
            Optimization result
        """
        cfg = self.config
        bounds = self.runner.get_active_bounds()

        # Create working directory
        work_dir = cfg.runs_root / f"direction_{direction_id:03d}"
        work_dir.mkdir(parents=True, exist_ok=True)

        # Initial point
        x0 = self._get_initial_point(direction_id, bounds)
        lb = np.array([b[0] for b in bounds])
        ub = np.array([b[1] for b in bounds])

        # History records
        cost_history: List[float] = []
        params_history: List[List[float]] = []

        # Save progress file
        progress_file = work_dir / "progress.json"

        def residual_func(params: np.ndarray) -> np.ndarray:
            """Residual function (Bayes unified normalized residuals)"""
            power_points = exp_data[:, 0]
            observations = exp_data[:, 1:4]

            predictions = self.runner.run(params, power_points)

            weights = np.asarray(cfg.output_weights, dtype=float)
            residuals = compute_bayes_normalized_residuals(
                predictions,
                observations,
                weights=weights,
            ).flatten()

            cost, _ = compute_bayes_nrmse_percent(
                predictions,
                observations,
                weights=weights,
            )
            cost_history.append(cost)
            params_history.append(params.tolist())

            # Update progress
            self._save_progress(
                progress_file,
                direction_id,
                len(cost_history),
                cost,
                min(cost_history),
            )

            if len(cost_history) % cfg.log_every_eval == 0:
                print(f"    [eval {len(cost_history)}] NRMSE = {cost:.4f}%")

            return residuals

        # Run optimization
        result = least_squares(
            residual_func,
            x0,
            bounds=(lb, ub),
            max_nfev=cfg.max_nfev,
            ftol=cfg.ftol,
            xtol=cfg.xtol,
            diff_step=cfg.diff_step,
        )

        final_predictions = self.runner.run(result.x, exp_data[:, 0])
        final_cost, _ = compute_bayes_nrmse_percent(
            final_predictions,
            exp_data[:, 1:4],
            weights=np.asarray(cfg.output_weights, dtype=float),
        )

        # Save results
        result_dict = {
            "direction_id": direction_id,
            "params": result.x.tolist(),
            "cost": final_cost,
            "nfev": result.nfev,
            "status": int(result.status),
            "message": result.message,
            "x0": x0.tolist(),
            "cost_history": cost_history,
        }

        result_file = work_dir / "result.json"
        with open(result_file, "w") as f:
            json.dump(result_dict, f, indent=2)

        return result_dict

    def _get_initial_point(
        self,
        direction_id: int,
        bounds: List,
    ) -> np.ndarray:
        """
        Generate initial point

        Parameters
        ----------
        direction_id : int
            Direction ID
        bounds : list
            Parameter bounds

        Returns
        -------
        np.ndarray
            Initial parameters
        """
        cfg = self.config
        n_params = len(bounds)

        if cfg.init_strategy == "legacy_single":
            # All directions use midpoint
            return np.array([(b[0] + b[1]) / 2 for b in bounds])

        elif cfg.init_strategy == "legacy_multi":
            # Random initialization
            rng = np.random.default_rng(seed=1000 + direction_id)
            return np.array([rng.uniform(b[0], b[1]) for b in bounds])

        elif cfg.init_strategy == "directional":
            # Start from midpoint, offset along a random direction
            center = np.array([(b[0] + b[1]) / 2 for b in bounds])
            scales = np.array([(b[1] - b[0]) / 2 for b in bounds])

            rng = np.random.default_rng(seed=1000 + direction_id)
            direction = rng.normal(0, 1, size=n_params)
            direction = direction / (np.linalg.norm(direction) + 1e-10)

            x0 = center + cfg.direction_scale * scales * direction

            # Clip to within bounds
            for i in range(n_params):
                x0[i] = np.clip(x0[i], bounds[i][0], bounds[i][1])

            return x0

        else:
            raise ValueError(f"Unknown initialization strategy: {cfg.init_strategy}")

    def _save_progress(
        self,
        filepath: Path,
        direction_id: int,
        n_eval: int,
        current_cost: float,
        best_cost: float,
    ) -> None:
        """Save progress information"""
        data = {
            "direction_id": direction_id,
            "n_eval": n_eval,
            "current_cost": current_cost,
            "best_cost": best_cost,
            "ts": time.time(),
        }
        with open(filepath, "w") as f:
            json.dump(data, f)

    def _monitor_loop(self, stop_event: threading.Event) -> None:
        """Monitoring thread"""
        cfg = self.config
        poll_interval = max(0.1, cfg.monitor_interval_sec)

        while not stop_event.is_set():
            # Collect progress from all directions
            rows = []
            for i in range(1, cfg.n_directions + 1):
                progress_file = cfg.runs_root / f"direction_{i:03d}" / "progress.json"
                if progress_file.exists():
                    try:
                        with open(progress_file) as f:
                            d = json.load(f)
                        rows.append((
                            i,
                            d.get("n_eval", 0),
                            d.get("current_cost"),
                            d.get("best_cost"),
                        ))
                    except Exception:
                        pass

            # Summary info can be printed here (currently omitted to reduce output)

            time.sleep(poll_interval)
