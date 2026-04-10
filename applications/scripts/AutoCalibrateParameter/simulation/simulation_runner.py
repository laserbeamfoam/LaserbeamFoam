"""
Simulation Runner

Wraps OpenFOAMCaseManager, providing a simplified interface for the optimization loop
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Callable, List, Optional, Tuple

import numpy as np

from ..models.foam_case_manager import OpenFOAMCaseManager
from ..utils.param_filter import filter_params_for_optimization, merge_active_and_fixed_params
from ..utils.param_registry import PARAM_NAMES, get_param_bounds_from_config
from .objective import compute_bayes_normalized_residuals, compute_bayes_nrmse_percent

if TYPE_CHECKING:
    from ..config.base_config import BaseConfig


class SimulationRunner:
    """
    Unified simulation runner

    Wraps OpenFOAMCaseManager, providing a simplified interface for the optimization loop

    Attributes
    ----------
    config : BaseConfig
        Configuration object
    case_manager : OpenFOAMCaseManager
        OpenFOAM case manager
    work_dir : Path
        Working directory
    eval_count : int
        Evaluation counter
    """

    def __init__(
        self,
        config: "BaseConfig",
        work_dir: Optional[Path] = None,
    ):
        """
        Initialize the simulation runner

        Parameters
        ----------
        config : BaseConfig
            Configuration object
        work_dir : Path, optional
            Working directory; defaults to config.runs_root / "simulations"
        """
        self.config = config
        self.work_dir = work_dir or config.runs_root / "simulations"
        self.work_dir = Path(self.work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

        # Initialize the OpenFOAM case manager
        self.case_manager = OpenFOAMCaseManager.from_config(config)

        # Parameter filtering configuration
        self.all_param_names = PARAM_NAMES
        config_bounds = get_param_bounds_from_config(config)

        self.active_names, self.active_bounds, self.fixed_dict = filter_params_for_optimization(
            PARAM_NAMES,
            config_bounds,
            config.active_params,
            config.fixed_values,
        )

        print(f"Optimization parameters: {self.active_names}")
        if self.fixed_dict:
            print(f"Fixed parameters: {self.fixed_dict}")

        # Evaluation counter
        self.eval_count = 0

    def get_active_bounds(self) -> List[Tuple[float, float]]:
        """Returns the bounds of the optimization parameters (for use by optimizers)"""
        return self.active_bounds

    def run(
        self,
        params: np.ndarray,
        power_points: np.ndarray,
        job_id: Optional[str | int] = None,
        on_power_complete: Optional[Callable[[int, float, np.ndarray], None]] = None,
        start_index: int = 0,
    ) -> np.ndarray:
        """
        Run simulation and return predictions

        Parameters
        ----------
        params : np.ndarray, shape (n_active,)
            Optimization parameter values (active parameters only)
        power_points : np.ndarray
            List of power points [P1, P2, ...]
        job_id : str | int, optional
            Job ID used for archiving results
        on_power_complete : Callable[[int, float, np.ndarray], None], optional
            Callback function called after each power point completes
            Arguments: (power_index, power_value, predictions_so_far)
        start_index : int, optional
            Starting power point index (for resume from checkpoint), defaults to 0

        Returns
        -------
        predictions : np.ndarray, shape (n_powers, 3)
            Predictions [width, depth, area] in units: μm, μm, μm²
        """
        params = np.asarray(params).flatten()
        power_points = np.asarray(power_points).flatten()

        # Merge active + fixed -> full parameter vector
        full_params = merge_active_and_fixed_params(
            self.all_param_names,
            self.active_names,
            params,
            self.fixed_dict,
        )

        sigma = full_params[self.all_param_names.index("sigma")]
        marangoni = full_params[self.all_param_names.index("marangoni")]
        substrate_temp = full_params[self.all_param_names.index("substrate_temp")]
        absorptivity = full_params[self.all_param_names.index("absorptivity")]
        recoil_coeff = full_params[self.all_param_names.index("recoilCoeff")]
        radius_flavour = full_params[self.all_param_names.index("radius_flavour")]
        laser_radius = full_params[self.all_param_names.index("laser_radius")]

        n_powers = len(power_points)
        predictions = np.zeros((n_powers, 3))

        for i, power in enumerate(power_points):
            self.eval_count += 1

            try:
                # Update parameters
                self.case_manager.update_parameters(
                    sigma,
                    marangoni,
                    substrate_temp,
                    absorptivity,
                    recoil_coeff=recoil_coeff,
                    radius_flavour=radius_flavour,
                    laser_radius=laser_radius,
                )
                self.case_manager.set_power(power, absorptivity)

                # Run simulation
                self.case_manager.run_simulation()
                self.case_manager.run_postprocess()

                # Read results (m -> μm)
                metrics = self.case_manager.read_metrics()
                predictions[i, 0] = metrics["width_mean_m"] * 1e6
                predictions[i, 1] = metrics["depth_mean_m"] * 1e6
                predictions[i, 2] = metrics["area_mean_m2"] * 1e12

                print(
                    f"  [eval {self.eval_count}] P={power:.0f}W: "
                    f"w={predictions[i, 0]:.2f}μm, d={predictions[i, 1]:.2f}μm"
                )

            except Exception as e:
                print(f"  [eval {self.eval_count}] P={power:.0f}W simulation failed: {e}")
                predictions[i, :] = np.nan

            finally:
                # Archive results
                if job_id is not None:
                    archive_dir = self.config.runs_root / f"{job_id}" / f"{int(power)}W"
                    archive_dir.mkdir(parents=True, exist_ok=True)
                    archive_mode = getattr(self.config, "archive_mode", "latest")
                    self.case_manager.archive_results(archive_dir, mode=archive_mode)

                # Call callback function (if provided)
                if on_power_complete is not None:
                    try:
                        global_idx = start_index + i
                        on_power_complete(global_idx, power, predictions[: i + 1].copy())
                    except Exception as e:
                        print(f"  [Warning] Power point completion callback failed: {e}")

        return predictions

    def run_single(
        self,
        params: np.ndarray,
        power: float,
    ) -> Tuple[float, float, float]:
        """
        Run simulation for a single power point

        Parameters
        ----------
        params : np.ndarray
            Parameters (active parameter vector)
        power : float
            Laser power (W)

        Returns
        -------
        tuple
            (width, depth, area) in units: μm, μm, μm²
        """
        result = self.run(params, np.array([power]))
        return tuple(result[0])

    def evaluate_objective(
        self,
        params: np.ndarray,
        exp_data: np.ndarray,
    ) -> float:
        """
        Evaluate the objective function (Bayes unified NRMSE%)

        Parameters
        ----------
        params : np.ndarray
            Parameters (active parameter vector)
        exp_data : np.ndarray, shape (n, 4)
            Experimental data [power, width, depth, area]

        Returns
        -------
        float
            NRMSE percentage
        """
        power_points = exp_data[:, 0]
        observations = exp_data[:, 1:4]

        predictions = self.run(params, power_points)
        cost, _ = compute_bayes_nrmse_percent(
            predictions,
            observations,
            weights=np.array(self.config.output_weights, dtype=float),
        )
        return cost

    def compute_residuals_flat(
        self,
        params: np.ndarray,
        exp_data: np.ndarray,
    ) -> np.ndarray:
        """
        Compute the flattened residual vector (for use with scipy.optimize.least_squares)

        Returns the normalized residuals corresponding to the Bayes unified objective
        (supports output_weights).
        """
        power_points = exp_data[:, 0]
        observations = exp_data[:, 1:4]

        predictions = self.run(params, power_points)

        residuals = compute_bayes_normalized_residuals(
            predictions,
            observations,
            weights=np.array(self.config.output_weights, dtype=float),
        )
        return residuals.flatten()

    def reset_eval_count(self) -> None:
        """Reset the evaluation counter"""
        self.eval_count = 0
