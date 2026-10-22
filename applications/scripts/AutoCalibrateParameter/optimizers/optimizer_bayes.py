"""
Bayesian optimizer

Uses scikit-optimize for Bayesian optimization (point estimate)
"""

from __future__ import annotations

import json
import os
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Tuple

import numpy as np

try:
    from skopt import Optimizer
    from skopt.space import Real
except ImportError:
    Optimizer = None
    Real = None

try:
    from pyDOE import lhs
except ImportError:
    lhs = None

import pandas as pd

from .base_optimizer import BaseOptimizer, OptimizationResult
from ..simulation.objective import compute_bayes_nrmse_percent
from ..utils.param_registry import (
    get_param_csv_candidates,
    get_param_csv_column,
    get_param_csv_column_map,
)

if TYPE_CHECKING:
    from ..config.bayes_config import BayesConfig
    from ..simulation.simulation_runner import SimulationRunner


def _evaluate_single_point(args: Tuple) -> Tuple[int, List[float], float]:
    """
    Single-point evaluation function (for parallel use)

    Parameters
    ----------
    args : tuple
        (job_id, params, config_dict, exp_data, output_weights)

    Returns
    -------
    tuple
        (job_id, params, cost)
    """
    job_id, params, config_dict, exp_data, output_weights = args

    # Rebuild configuration and runner
    from ..config.bayes_config import BayesConfig
    from ..simulation.simulation_runner import SimulationRunner
    from ..simulation.objective import compute_bayes_nrmse_percent

    config = BayesConfig(**config_dict)
    runner = SimulationRunner(config)

    # Evaluate
    power_points = exp_data[:, 0]
    observations = exp_data[:, 1:4]

    predictions = runner.run(np.array(params), power_points, job_id=job_id)

    cost, _ = compute_bayes_nrmse_percent(
        predictions,
        observations,
        weights=np.asarray(output_weights, dtype=float),
    )

    return job_id, params, cost


class BayesianOptimizer(BaseOptimizer):
    """
    Bayesian optimizer (scikit-optimize)

    Features:
    - Gaussian process surrogate + acquisition function
    - Returns point estimate (optimal solution)
    - Supports batch parallel evaluation

    Attributes
    ----------
    config : BayesConfig
        Bayesian optimization configuration
    runner : SimulationRunner
        Simulation runner
    """

    def __init__(self, config: "BayesConfig", runner: "SimulationRunner"):
        """
        Initialize the Bayesian optimizer

        Parameters
        ----------
        config : BayesConfig
            Bayesian optimization configuration
        runner : SimulationRunner
            Simulation runner
        """
        if Optimizer is None:
            raise ImportError("scikit-optimize is required: pip install scikit-optimize")

        super().__init__(config)
        self.config = config
        self.runner = runner

        # Output weights
        self.output_weights = np.array(config.output_weights)

    def get_name(self) -> str:
        return "Bayesian"

    def _get_history_csv_path(self) -> Path:
        """Get the unique history CSV path (always stored under runs_root)."""
        history_name = "bayes_history.csv"
        if getattr(self.config, "warmstart_csv", None):
            configured_raw = str(self.config.warmstart_csv)
            configured_name = Path(configured_raw).name
            if configured_raw != configured_name:
                print(f"Warning: warmstart_csv='{configured_raw}' contains a path; normalized to filename '{configured_name}' stored in runs_root")
            if configured_name:
                history_name = configured_name
        return self.config.runs_root / history_name

    def _get_progress_file_path(self, job_id: int) -> Path:
        """Get the progress file path for the specified job (always stored under runs_root)."""
        return self.config.runs_root / f"bayes_progress_job{job_id}.csv"

    def _migrate_legacy_files(self) -> None:
        """Migrate legacy case root directory files to runs_root, keeping only one copy under runs."""
        self.config.runs_root.mkdir(parents=True, exist_ok=True)

        # Migrate bayes_history.csv (if it does not exist in runs)
        history_in_runs = self._get_history_csv_path()
        legacy_history = self.config.case_dir / history_in_runs.name
        if legacy_history.exists() and legacy_history != history_in_runs:
            if history_in_runs.exists():
                try:
                    legacy_history.unlink()
                    print(f"Cleaned up legacy history file: {legacy_history}")
                except Exception as e:
                    print(f"Warning: legacy history {legacy_history} detected but cannot be deleted ({e}); continuing with {history_in_runs}")
            else:
                shutil.move(str(legacy_history), str(history_in_runs))
                print(f"Migrated history file to runs: {history_in_runs}")

        # Migrate bayes_progress_job*.csv
        for legacy_progress in sorted(self.config.case_dir.glob("bayes_progress_job*.csv")):
            target_progress = self.config.runs_root / legacy_progress.name
            if legacy_progress == target_progress:
                continue
            if target_progress.exists():
                try:
                    legacy_progress.unlink()
                    print(f"Cleaned up legacy progress file: {legacy_progress}")
                except Exception as e:
                    print(f"Warning: cannot clean up legacy progress file {legacy_progress}: {e}")
            else:
                shutil.move(str(legacy_progress), str(target_progress))
                print(f"Migrated progress file to runs: {target_progress}")

    def optimize(self, exp_data: np.ndarray) -> OptimizationResult:
        """
        Execute Bayesian optimization

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
        print(f"Bayesian optimization (scikit-optimize)")
        print(f"{'='*60}")
        print(f"Acquisition function: {cfg.acq_func}")
        print(f"Initial points: {cfg.n_initial_points}, batches: {cfg.n_batches}, batch size: {cfg.batch_size}")
        print()

        # Define search space (using the active parameter bounds from runner)
        bounds = self.runner.get_active_bounds()
        param_names = self.runner.active_names
        dimensions = [
            Real(low, high, name=name)
            for (low, high), name in zip(bounds, param_names)
        ]

        # Initialize the skopt optimizer
        optimizer = Optimizer(
            dimensions=dimensions,
            base_estimator="GP",
            acq_func=cfg.acq_func,
            acq_func_kwargs={"xi": cfg.xi, "kappa": cfg.kappa},
            n_initial_points=0,  # We handle initial sampling ourselves
            random_state=42,
        )

        all_x: List[List[float]] = []
        all_y: List[float] = []
        all_predictions: List[Optional[np.ndarray]] = []  # Store all prediction results
        best_cost = float("inf")

        # Get power points for saving CSV
        power_points = exp_data[:, 0]

        # Migrate all legacy files to runs_root
        self._migrate_legacy_files()

        # Check for warm start from CSV under runs
        history_csv_path = self._get_history_csv_path()
        if history_csv_path.exists():
            try:
                # Load history and predictions
                warmstart_x, warmstart_y, warmstart_preds = self._load_history_csv(str(history_csv_path))
                if warmstart_x and warmstart_y:
                    all_x.extend(warmstart_x)
                    all_y.extend(warmstart_y)
                    # Use loaded predictions instead of None
                    all_predictions.extend(warmstart_preds)
                    best_cost = min(all_y)
                    print(f"Warm start from CSV: {len(warmstart_y)} samples, best RMSE = {best_cost:.2f}")
            except Exception as e:
                print(f"Failed to load CSV warm start data: {e}")

        # Check if resuming from history.json (Deprecated - use CSV only)
        # history_file = cfg.runs_root / "history.json"
        # if cfg.resume and history_file.exists():
        #     try:
        #         with open(history_file, "r") as f:
        #             history_data = json.load(f)
        #         resume_x = history_data.get("x_iters", [])
        #         resume_y = history_data.get("y_iters", [])
        #         if resume_x and resume_y:
        #             # Merge data (avoid duplicates)
        #             existing_set = set(tuple(x) for x in all_x)
        #             for x, y in zip(resume_x, resume_y):
        #                 if tuple(x) not in existing_set:
        #                     all_x.append(x)
        #                     all_y.append(y)
        #                     all_predictions.append(None)  # No predictions for resumed data
        #                     existing_set.add(tuple(x))
        #             if all_y:
        #                 best_cost = min(all_y)
        #             print(f"Resumed from history.json: {len(all_y)} samples, best RMSE = {best_cost:.2f}")
        #     except Exception as e:
        #         print(f"Failed to resume history: {e}")

        # Provide historical data to the optimizer
        if all_x and all_y:
            # Important: Ensure NO NaNs in all_y before telling optimizer
            valid_indices = [i for i, y in enumerate(all_y) if not np.isnan(y)]
            if len(valid_indices) < len(all_y):
                 print(f"Warning: ignoring {len(all_y) - len(valid_indices)} records with NaN objective values")
            
            clean_x = [all_x[i] for i in valid_indices]
            clean_y = [all_y[i] for i in valid_indices]

            # Clip historical points to current bounds (handles cases where bounds were widened/narrowed)
            bounds_list = [space_dim.bounds for space_dim in optimizer.space.dimensions]
            clipped_x = []
            n_clipped = 0
            for pt in clean_x:
                clipped_pt = [float(np.clip(v, lo, hi)) for v, (lo, hi) in zip(pt, bounds_list)]
                if clipped_pt != list(pt):
                    n_clipped += 1
                clipped_x.append(clipped_pt)
            if n_clipped > 0:
                print(f"Warning: {n_clipped} historical points exceeded new bounds and were clipped")

            if clipped_x and clean_y:
                optimizer.tell(clipped_x, clean_y)
                print(f"Total historical data: {len(clean_y)} valid samples")

        # Phase 1: Initial sampling
        n_initial_needed = max(0, cfg.n_initial_points - len(all_y))

        # Detect incomplete progress files and resume them first
        incomplete_jobs = self._detect_incomplete_jobs()
        if incomplete_jobs:
            print(f"\nDetected {len(incomplete_jobs)} incomplete jobs; resuming them first")
            for job_info in incomplete_jobs:
                job_id = job_info['job_id']
                params = job_info['params']
                print(f"  [Job {job_id}] Resuming evaluation: {params}")

                # Directly evaluate this incomplete job
                n_before = len(all_y)
                batch_results = self._evaluate_batch(
                    [params], exp_data, n_before,
                    all_x=all_x, all_y=all_y, all_predictions=all_predictions,
                    power_points=exp_data[:, 0],
                    explicit_job_ids=[job_id]
                )

                # Update history
                for param_vals, cost, pred in batch_results:
                    all_x.append(param_vals)
                    all_y.append(cost)
                    all_predictions.append(pred)
                    if not np.isnan(cost) and (best_cost is None or cost < best_cost):
                        best_cost = cost

                # Reduce the number of samples needed
                n_initial_needed = max(0, n_initial_needed - 1)

        if n_initial_needed > 0:
            print(f"\nPhase 1: Initial sampling ({n_initial_needed} points)")
            initial_samples = self._initial_sampling(n_initial_needed, bounds)

            for batch_start in range(0, n_initial_needed, cfg.batch_size):
                batch_end = min(batch_start + cfg.batch_size, n_initial_needed)
                batch_samples = initial_samples[batch_start:batch_end]

                # Record count before batch starts, to capture newly added results
                n_before = len(all_y)

                # Call _evaluate_batch; it saves each simulation result in real time
                batch_results = self._evaluate_batch(
                    batch_samples, exp_data, n_before,
                    all_x=all_x, all_y=all_y, all_predictions=all_predictions,
                    power_points=power_points
                )

                # Update best value and log
                for i, (params, cost, preds) in enumerate(batch_results):
                    if cost < best_cost:
                        best_cost = cost
                    self._log_iteration(np.array(params), cost, n_before + i + 1)

                print(f"  Batch complete, current best RMSE = {best_cost:.2f}")
                self._generate_interim_plots(
                    all_x, all_y, all_predictions, exp_data, power_points,
                    n_initial=cfg.n_initial_points,
                )

        # Phase 2: Bayesian optimization
        print(f"\nPhase 2: Bayesian optimization ({cfg.n_batches} batches)")

        # Re-tell only valid data just in case
        # (Though we did it before, newly added initial samples are automatically handled by loop or manual tell?)
        # skopt optimizer handles tell() incrementally. We don't need to re-tell everything.
        # But we DO need to tell the INITIAL SAMPLES we just ran if they weren't told.
        # However, we only told the `warmstart` data. The new `initial_samples` need to be told.
        # Since we use `ask()` in loop below, we must be careful.
        # Actually, standard usage is: tell initial points, then ask-tell loop.
        
        # Tell the initial samples we just ran
        if n_initial_needed > 0:
             # Identify the indices of the new samples in all_x/all_y
             # They are at the end: from len(all_y)-n_initial_needed to end?
             # No, we might have had failures.
             # Simply: we haven't told the optimizer about the new samples yet.
             # In skopt, `tell` returns a result object, but we manage state ourselves too.
             
             # We need to tell the optimizer about the samples obtained in Phase 1
             # We can just tell everything again? No, that duplicates.
             # We should tell the NEW samples.
             
             # The warmstart data was told. Now we tell the Phase 1 data.
             # Phase 1 added samples to all_x, all_y.
             # So we tell the slice that wasn't told.
             start_idx = len(all_x) - n_initial_needed # Assuming all succeeded?
             # That assumption is risky. Let's rely on checking what's new.
             # Actually, simpler: construct a clean list of EVERYTHING and re-initialize optimizer? 
             # No, optimizer has internal state (GP fit).
             
             # Better: just tell the new points from Phase 1.
             # Since we know exactly which ones they are (from batch_results).
             # We didn't accumulate batch_results from the loop easily.
             
             # Retrospective fix: We iterate batches in Phase 1. 
             # We should collect them and tell.
             pass 
             
             # Actually, let's just make sure we tell ALL valid points that haven't been told.
             # But optimizer.tell() is stateful.
             # Let's clean this up:
             # 1. New Optimizer is created.
             # 2. We tell warmstart data.
             # 3. We run Phase 1 (Initial Sampling).
             # 4. We MUST tell Phase 1 data to optimizer before Phase 2.
        
        # Let's grab the data generated in Phase 1
        # It's all_x/all_y minus the warmstart count.
        # But wait, we don't know the exact count easily unless we tracked it.
        # Let's count how many we told initially.
        n_told = 0
        if all_x and all_y:
             # We told everything available at start
             # But we filtered NaNs.
             pass 
        
        # To be safe and simple: 
        # For the Phase 1 points, we should iterate and tell them.
        # But doing it inside the loop is messy.
        
        # Correct approach:
        # We collected all_x/all_y. The optimizer only knows about the warmstart part.
        # We should tell the rest.
        
        # Let's re-fit the optimizer with ALL current data to be sure.
        # But skopt doesn't have a 'reset' easily.
        # Actually, we can just tell the new points.
        
        # IMPORTANT: The original code didn't tell Phase 1 points to optimizer explicitly!
        # It relied on `optimizer` being used in Phase 2? No.
        # If we don't tell Phase 1 points, Phase 2 starts with only warmstart info.
        # That's a bug in original code too, potentially.
        # Wait, usually `n_initial_points` in `Optimizer` handles this?
        # We set `n_initial_points=0` and do it manually.
        
        # So we MUST tell the points from Phase 1.
        # Let's identify them.
        # We can just check `optimizer.Xi` to see what it has (if accessible) or just track manually.
        
        # Simplest hack: We know `all_x` contains everything.
        # We can create a NEW optimizer for Phase 2, fed with ALL history.
        optimizer = Optimizer(
            dimensions=dimensions,
            base_estimator="GP",
            acq_func=cfg.acq_func,
            acq_func_kwargs={"xi": cfg.xi, "kappa": cfg.kappa},
            n_initial_points=0,
            random_state=42,
        )
        
        # Tell ALL valid history
        valid_indices = [i for i, y in enumerate(all_y) if not np.isnan(y)]
        clean_x = [all_x[i] for i in valid_indices]
        clean_y = [all_y[i] for i in valid_indices]

        # Clip to current bounds (handles cases where bounds were widened/narrowed)
        bounds_list = [d.bounds for d in optimizer.space.dimensions]
        clean_x = [[float(np.clip(v, lo, hi)) for v, (lo, hi) in zip(pt, bounds_list)] for pt in clean_x]

        if clean_x and clean_y:
            optimizer.tell(clean_x, clean_y)


        for batch_idx in range(cfg.n_batches):
            # Get candidate points
            candidates = optimizer.ask(n_points=cfg.batch_size)
            candidates = [list(c) for c in candidates]

            print(f"\nBatch {batch_idx + 1}/{cfg.n_batches}")

            # Record count before batch starts
            n_before = len(all_y)

            # Call _evaluate_batch; it saves each simulation result in real time
            batch_results = self._evaluate_batch(
                candidates, exp_data, n_before,
                all_x=all_x, all_y=all_y, all_predictions=all_predictions,
                power_points=power_points
            )

            # Collect batch results to inform the optimizer
            batch_x = [r[0] for r in batch_results]
            batch_y = [r[1] for r in batch_results]

            # Update best value and log
            for i, (params, cost, preds) in enumerate(batch_results):
                if cost < best_cost:
                    best_cost = cost
                self._log_iteration(np.array(params), cost, n_before + i + 1)

            # Only tell valid results
            valid_batch_indices = [i for i, y in enumerate(batch_y) if not np.isnan(y)]
            if valid_batch_indices:
                clean_batch_x = [batch_x[i] for i in valid_batch_indices]
                clean_batch_y = [batch_y[i] for i in valid_batch_indices]
                optimizer.tell(clean_batch_x, clean_batch_y)

            if len(valid_batch_indices) < len(batch_y):
                 print(f"  Warning: ignoring {len(batch_y) - len(valid_batch_indices)} failed/NaN results in this batch")

            print(f"  Batch {batch_idx + 1} complete, best RMSE = {best_cost:.2f}")
            self._generate_interim_plots(
                all_x, all_y, all_predictions, exp_data, power_points,
                n_initial=cfg.n_initial_points,
            )

        # Extract the optimal solution
        best_idx = int(np.argmin(all_y))
        best_params = np.array(all_x[best_idx])

        return OptimizationResult(
            method="Bayesian",
            best_params=best_params,
            best_cost=all_y[best_idx],
            n_evaluations=len(all_y),
            history={"params": all_x, "costs": all_y},
            message=f"Bayesian optimization complete, {len(all_y)} evaluations total",
        )

    def _generate_interim_plots(
        self,
        all_x: List[List[float]],
        all_y: List[float],
        all_predictions: List,
        exp_data: np.ndarray,
        power_points: np.ndarray,
        n_initial: int = 0,
    ) -> None:
        """Generate intermediate visualization plots after each batch (overwrite update, no file accumulation)"""
        from ..postprocess.comparison import save_comparison_plot, plot_predictions_vs_experiments
        from ..postprocess.convergence import save_convergence_plot

        valid_indices = [i for i, y in enumerate(all_y) if not np.isnan(y)]
        if not valid_indices:
            return

        best_idx = min(valid_indices, key=lambda i: all_y[i])
        best_cost = all_y[best_idx]
        best_preds = all_predictions[best_idx]

        if best_preds is None:
            return

        # Take only width/depth columns, ignore area
        observations = exp_data[:, 1:3]
        best_preds_wd = best_preds[:, :2]

        plots_dir = self.config.runs_root / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)

        n_evals = len(valid_indices)

        try:
            save_comparison_plot(
                plots_dir,
                power_points,
                observations,
                best_preds_wd,
                output_names=["Width", "Depth"],
                output_units=["um", "um"],
                nrmse=best_cost,
                title=f"Best result (n={n_evals}, NRMSE={best_cost:.2f}%)",
            )
            plot_predictions_vs_experiments(
                plots_dir,
                observations,
                best_preds_wd,
                output_names=["Width", "Depth"],
                output_units=["um", "um"],
                title=f"Predictions vs Experiments (n={n_evals})",
            )
        except Exception as e:
            print(f"  [Warning] Failed to generate comparison plot: {e}")

        try:
            valid_costs = [all_y[i] for i in valid_indices]
            save_convergence_plot(
                plots_dir, valid_costs, n_initial=n_initial,
                title="Optimization Convergence",
            )
        except Exception as e:
            print(f"  [Warning] Failed to generate convergence plot: {e}")

        try:
            from ..postprocess.convergence import plot_parameter_evolution
            valid_params = [all_x[i] for i in valid_indices]
            plot_parameter_evolution(
                plots_dir,
                [np.array(p) for p in valid_params],
                self.runner.active_names,
                title="Parameter Evolution",
            )
        except Exception as e:
            print(f"  [Warning] Failed to generate parameter evolution plot: {e}")

        print(f"  [Plots] Updated: {plots_dir}")

    def _initial_sampling(
        self,
        n_points: int,
        bounds: List[Tuple[float, float]],
    ) -> List[List[float]]:
        """Initial LHS sampling"""
        rng_seed = 42 if self.config.seed is None else int(self.config.seed)

        if lhs is None:
            # Fallback: random sampling
            np.random.seed(rng_seed)
            samples = []
            for _ in range(n_points):
                sample = [
                    np.random.uniform(low, high)
                    for low, high in bounds
                ]
                samples.append(sample)
            return samples

        np.random.seed(rng_seed)
        lhs_samples = lhs(len(bounds), n_points)
        samples = []
        for row in lhs_samples:
            sample = [
                bounds[i][0] + row[i] * (bounds[i][1] - bounds[i][0])
                for i in range(len(bounds))
            ]
            samples.append(sample)
        return samples

    def _evaluate_batch(
        self,
        candidates: List[List[float]],
        exp_data: np.ndarray,
        offset: int,
        all_x: List[List[float]] = None,
        all_y: List[float] = None,
        all_predictions: List = None,
        power_points: np.ndarray = None,
        explicit_job_ids: List[int] = None,
    ) -> List[Tuple[List[float], float, Optional[np.ndarray]]]:
        """
        Evaluate a batch of candidate points
        """
        results = []

        for idx, params in enumerate(candidates):
            if explicit_job_ids:
                job_id = explicit_job_ids[idx]
            else:
                job_id = offset + idx + 1

            # Check for a partially completed progress file
            progress_file = self._get_progress_file_path(job_id)
            partial_preds = None
            start_power_idx = 0

            if progress_file.exists():
                try:
                    partial_preds, start_power_idx = self._load_partial_progress(
                        job_id, params, power_points
                    )
                    if partial_preds is not None and start_power_idx > 0:
                        print(f"  [Job {job_id}] Partial progress detected; resuming from power point {start_power_idx+1}/{len(power_points)}")
                except Exception as e:
                    print(f"  [Job {job_id}] Cannot load partial progress: {e}; starting from scratch")
                    partial_preds = None
                    start_power_idx = 0
            # Dynamically build parameter log string
            param_strs = []
            for i, name in enumerate(self.runner.active_names):
                if name == "sigma":
                    param_strs.append(f"sigma={params[i]:.4f}")
                elif name == "marangoni":
                    param_strs.append(f"marangoni={params[i]:.2e}")
                elif name == "substrate_temp":
                    param_strs.append(f"T_s={params[i]:.1f}K")
                elif name == "absorptivity":
                    param_strs.append(f"absorptivity={params[i]:.3f}")
                elif name == "recoilCoeff":
                    param_strs.append(f"recoilCoeff={params[i]:.3f}")
                elif name == "radius_flavour":
                    param_strs.append(f"Radius_Flavour={params[i]:.3f}")
                else:
                    param_strs.append(f"{name}={params[i]:.4g}")
            print(f"  [Job {job_id}] Parameters: {', '.join(param_strs)}")

            try:
                exp_power_points = exp_data[:, 0]
                observations = exp_data[:, 1:4]

                # Create callback function to save intermediate results after each power point
                partial_predictions_storage = {'data': None}

                def on_power_complete(power_idx: int, power_val: float, preds_so_far: np.ndarray):
                    """Callback invoked after each power point completes"""
                    partial_predictions_storage['data'] = preds_so_far.copy()

                    # Compute partial NRMSE for completed power points (consistent with Bayes objective)
                    n_completed = power_idx + 1
                    partial_obs = observations[:n_completed, :]
                    partial_preds = preds_so_far

                    partial_cost, _ = compute_bayes_nrmse_percent(
                        partial_preds,
                        partial_obs,
                        weights=np.asarray(self.output_weights, dtype=float),
                    )

                    print(f"  [Job {job_id}] Power point {power_idx+1}/{len(exp_power_points)} complete "
                          f"(P={power_val:.0f}W), partial NRMSE={partial_cost:.2f}%")

                    # Save intermediate results to temporary CSV in real time
                    if all_x is not None and power_points is not None:
                        self._save_partial_progress(
                            job_id, params, partial_preds, power_points[:n_completed],
                            partial_cost, n_completed, len(exp_power_points)
                        )

                # If partially completed data is available, resume from checkpoint
                if partial_preds is not None and start_power_idx > 0:
                    # Initialize predictions array with already completed data
                    predictions = partial_preds.copy()

                    # Only run the remaining power points
                    remaining_powers = exp_power_points[start_power_idx:]
                    if len(remaining_powers) > 0:
                        # Create enhanced callback that passes the full prediction array
                        def enhanced_callback(power_idx: int, power_val: float, partial_remaining: np.ndarray):
                            """Wrapper callback that passes complete data"""
                            # Update the full predictions array
                            predictions[power_idx, :] = partial_remaining[-1, :]

                            # Call the original callback with complete data
                            if on_power_complete is not None:
                                on_power_complete(power_idx, power_val, predictions[:power_idx+1].copy())

                        remaining_preds = self.runner.run(
                            np.array(params), remaining_powers, job_id=job_id,
                            on_power_complete=enhanced_callback if on_power_complete else None,
                            start_index=start_power_idx  # Pass starting index
                        )
                        # Fill in the remaining predictions
                        predictions[start_power_idx:] = remaining_preds
                else:
                    # Run all power points from scratch
                    predictions = self.runner.run(
                        np.array(params), exp_power_points, job_id=job_id,
                        on_power_complete=on_power_complete
                    )

                cost, nrmse_per_output = compute_bayes_nrmse_percent(
                    predictions,
                    observations,
                    weights=np.asarray(self.output_weights, dtype=float),
                )
                output_names = ['width', 'depth', 'area']
                weights = np.asarray(self.output_weights, dtype=float)
                nrmse_str = ', '.join([
                    f"{name}={nrmse_per_output[i]:.2f}%{'*' if weights[i] > 0 else ''}"
                    for i, name in enumerate(output_names)
                ])
                print(f"  [Job {job_id}] NRMSE = {cost:.2f}% ({nrmse_str})")
                results.append((params, cost, predictions))

            except Exception as e:
                print(f"  [Job {job_id}] Failed: {e}")
                results.append((params, np.nan, None)) # Use np.nan for failure instead of 1e16 to handle properly

            # Save in real time
            if all_x is not None and all_y is not None and all_predictions is not None:
                all_x.append(results[-1][0])
                all_y.append(results[-1][1])
                all_predictions.append(results[-1][2])
                self._save_history(all_x, all_y)
                if power_points is not None:
                    self._save_detailed_csv(all_x, all_y, all_predictions, power_points)
                print(f"  [Job {job_id}] Data saved ({len(all_y)} records total)")

                # Delete temporary progress file (if it exists)
                progress_file = self._get_progress_file_path(job_id)
                if progress_file.exists():
                    try:
                        progress_file.unlink()
                        print(f"  [Job {job_id}] Temporary progress file deleted")
                    except Exception as e:
                        print(f"  [Warning] Cannot delete temporary progress file: {e}")

        return results

    def _save_history(
        self,
        all_x: List[List[float]],
        all_y: List[float],
    ) -> None:
        """Save optimization history (deprecated; no longer saves json)"""
        pass
        # history_file = self.config.runs_root / "history.json"
        # ... (removed)

    def _save_detailed_csv(
        self,
        all_x: List[List[float]],
        all_y: List[float],
        all_predictions: List[Optional[np.ndarray]],
        power_points: np.ndarray,
    ) -> None:
        """
        Save detailed optimization history to a CSV file (single file under runs_root only)
        """
        rows = []
        for i, (params, cost, preds) in enumerate(zip(all_x, all_y, all_predictions)):
            # Dynamically build parameter columns (based on active_names)
            row = {}
            for j, name in enumerate(self.runner.active_names):
                row[get_param_csv_column(name)] = params[j]

            # Add fixed parameters (to record the complete configuration)
            for name, value in self.runner.fixed_dict.items():
                row[get_param_csv_column(name)] = value

            if preds is not None:
                for j, power in enumerate(power_points):
                    row[f"width_um_{int(power)}W"] = preds[j, 0]
                    row[f"depth_um_{int(power)}W"] = preds[j, 1]
                    row[f"area_um2_{int(power)}W"] = preds[j, 2]
            else:
                for power in power_points:
                    row[f"width_um_{int(power)}W"] = np.nan
                    row[f"depth_um_{int(power)}W"] = np.nan
                    row[f"area_um2_{int(power)}W"] = np.nan

            row["objective"] = cost 
            row["status"] = "ok" if (cost is not None and not np.isnan(cost) and cost < 1e6) else "failed"

            rows.append(row)

        df = pd.DataFrame(rows)

        csv_file = self._get_history_csv_path()
        csv_file.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(csv_file, index=False)

    @staticmethod
    def _resolve_param_csv_column(df: pd.DataFrame, param_name: str) -> str:
        """Parse parameter column names (prefer new column names, compatible with legacy column names)."""
        for candidate in get_param_csv_candidates(param_name):
            if candidate in df.columns:
                return candidate
        raise ValueError(f"CSV missing parameter column: {param_name}, candidate columns={get_param_csv_candidates(param_name)}")

    def _save_partial_progress(
        self,
        job_id: int,
        params: List[float],
        predictions: np.ndarray,
        power_points: np.ndarray,
        partial_cost: float,
        n_completed: int,
        n_total: int,
    ) -> None:
        """
        Save partially completed progress to a temporary file

        Parameters
        ----------
        job_id : int
            Job ID
        params : List[float]
            Parameter list
        predictions : np.ndarray
            Completed prediction results
        power_points : np.ndarray
            Completed power points
        partial_cost : float
            Partial NRMSE
        n_completed : int
            Number of completed power points
        n_total : int
            Total number of power points
        """
        # Save to temporary progress file
        progress_file = self._get_progress_file_path(job_id)

        # Dynamically build parameter row (consistent with _save_detailed_csv)
        row = {"job_id": job_id}
        for j, name in enumerate(self.runner.active_names):
            row[get_param_csv_column(name)] = params[j]

        # Add fixed parameters
        for name, value in self.runner.fixed_dict.items():
            row[get_param_csv_column(name)] = value

        row["progress"] = f"{n_completed}/{n_total}"
        row["partial_nrmse"] = partial_cost
        row["timestamp"] = pd.Timestamp.now().strftime("%Y/%m/%d %H:%M")

        # Add completed power point data
        for i, power in enumerate(power_points):
            row[f"width_um_{int(power)}W"] = predictions[i, 0]
            row[f"depth_um_{int(power)}W"] = predictions[i, 1]
            row[f"area_um2_{int(power)}W"] = predictions[i, 2]

        df = pd.DataFrame([row])
        df.to_csv(progress_file, index=False)

        print(f"  [Job {job_id}] Progress saved to {progress_file.name}")

    def _detect_incomplete_jobs(self) -> List[Dict]:
        """
        Detect incomplete jobs (those with bayes_progress_job*.csv but not in bayes_history.csv)

        Returns
        -------
        list of dict
            Each dictionary contains {'job_id': int, 'params': list, 'progress': str}
        """
        incomplete_jobs = []

        # Load experimental data to get the total number of power points
        from ..data import load_experiment_data
        try:
            exp_data = load_experiment_data(self.config.exp_csv)
            total_power_points = len(exp_data)
            power_points = exp_data[:, 0]
        except Exception as e:
            print(f"  [Warning] Cannot load experimental data: {e}")
            return incomplete_jobs

        # Scan progress files under runs_root
        self.config.runs_root.mkdir(parents=True, exist_ok=True)
        for progress_file in sorted(self.config.runs_root.glob("bayes_progress_job*.csv")):
            try:
                df = pd.read_csv(progress_file)
                if len(df) == 0:
                    continue

                row = df.iloc[0]
                job_id = int(row['job_id'])

                # Extract parameters (in active_names order)
                params = []
                for param_name in self.runner.active_names:
                    csv_col = self._resolve_param_csv_column(df, param_name)
                    if csv_col in row:
                        params.append(float(row[csv_col]))
                    else:
                        raise ValueError(f"Progress file missing parameter: {param_name}")

                # Determine progress by checking actual power point data (not relying on progress field)
                completed = 0
                for i, power in enumerate(power_points):
                    w_col = f"width_um_{int(power)}W"
                    d_col = f"depth_um_{int(power)}W"
                    a_col = f"area_um2_{int(power)}W"

                    # Check whether this power point's data exists and is valid
                    if all(col in row for col in [w_col, d_col, a_col]):
                        if not any(pd.isna(row[col]) for col in [w_col, d_col, a_col]):
                            completed = i + 1
                        else:
                            break  # Stop at first invalid data entry
                    else:
                        break  # Stop at first missing column

                # Only resume partially completed jobs (0 < completed < total_power_points)
                if 0 < completed < total_power_points:
                    print(f"  [Incomplete job found] Job {job_id}: {completed}/{total_power_points} power points completed")
                    incomplete_jobs.append({
                        'job_id': job_id,
                        'params': params,
                        'progress': f"{completed}/{total_power_points}",
                        'progress_file': progress_file
                    })

            except Exception as e:
                print(f"  [Warning] Failed to read progress file {progress_file.name}: {e}")
                continue

        return incomplete_jobs

    def _load_partial_progress(
        self,
        job_id: int,
        params: List[float],
        power_points: np.ndarray,
    ) -> Tuple[Optional[np.ndarray], int]:
        """
        Load partially completed progress

        Parameters
        ----------
        job_id : int
            Job ID
        params : List[float]
            Current parameters (for validation)
        power_points : np.ndarray
            All power points

        Returns
        -------
        tuple
            (partial_predictions, start_index)
            - partial_predictions: completed prediction data, shape (n_powers, 3)
            - start_index: index of the next power point to run
        """
        import pandas as pd

        progress_file = self._get_progress_file_path(job_id)
        if not progress_file.exists():
            return None, 0

        df = pd.read_csv(progress_file)
        if len(df) == 0:
            return None, 0

        row = df.iloc[0]

        # Validate that parameters match (optional, to prevent resuming the wrong job)
        param_match = True
        for i, name in enumerate(self.runner.active_names):
            csv_name = self._resolve_param_csv_column(df, name)
            if csv_name in row and abs(row[csv_name] - params[i]) > 1e-6:
                param_match = False
                break

        if not param_match:
            print(f"  [Warning] Job {job_id} parameter mismatch; ignoring progress file")
            return None, 0

        # Parse progress field (format: "1/3" means 1 completed out of 3 total)
        progress_str = str(row.get("progress", "0/0"))
        try:
            # Try standard format "2/3"
            completed, total = map(int, progress_str.split("/"))
        except:
            # Try handling Excel date format "3-Jan" -> 1/3 (January 3 -> month=1, day=3)
            try:
                import re
                # Match "number-month" format
                match = re.match(r'(\d+)-([A-Za-z]+)', progress_str)
                if match:
                    day = int(match.group(1))  # day (total)
                    month_str = match.group(2).capitalize()  # month name

                    # Month mapping: Excel converts 1/3 to January 3 (3-Jan)
                    month_map = {
                        "Jan": 1, "Feb": 2, "Mar": 3, "Apr": 4, "May": 5, "Jun": 6,
                        "Jul": 7, "Aug": 8, "Sep": 9, "Oct": 10, "Nov": 11, "Dec": 12
                    }

                    if month_str in month_map:
                        completed = month_map[month_str]  # month = completed count
                        total = day  # day = total count
                        print(f"  [Info] Excel date format '{progress_str}' parsed as {completed}/{total}")
                    else:
                        raise ValueError(f"Unrecognized month: {month_str}")
                else:
                    raise ValueError(f"Cannot parse: {progress_str}")
            except Exception as e:
                print(f"  [Warning] Cannot parse progress field: {progress_str} ({e})")
                return None, 0

        if completed == 0 or completed >= len(power_points):
            return None, 0

        # Read completed prediction data
        predictions = np.zeros((len(power_points), 3))
        predictions[:] = np.nan  # Initialize to NaN

        n_loaded = 0
        for i, power in enumerate(power_points[:completed]):
            w_col = f"width_um_{int(power)}W"
            d_col = f"depth_um_{int(power)}W"
            a_col = f"area_um2_{int(power)}W"

            if w_col in row and d_col in row and a_col in row:
                predictions[i, 0] = row[w_col]
                predictions[i, 1] = row[d_col]
                predictions[i, 2] = row[a_col]
                n_loaded += 1

        if n_loaded != completed:
            print(f"  [Warning] Expected to load {completed} power points, actually loaded {n_loaded}")

        return predictions, completed

    def _load_history_csv(
        self,
        csv_path: str,
    ) -> Tuple[List[List[float]], List[float], List[Optional[np.ndarray]]]:
        """
        Load historical data from a CSV file

        Returns
        -------
        tuple
            (all_x, all_y, all_predictions)
        """
        df = pd.read_csv(csv_path)

        all_x = []
        all_y = []
        all_predictions = []

        csv_col_map = get_param_csv_column_map()

        # Check that required columns exist
        required_csv_cols = []
        for param_name in self.runner.active_names:
            if param_name not in csv_col_map:
                raise ValueError(f"Unknown parameter name: {param_name}")
            csv_col = self._resolve_param_csv_column(df, param_name)
            required_csv_cols.append(csv_col)

        if not all(col in df.columns for col in required_csv_cols):
            raise ValueError(f"CSV missing required parameter columns: {required_csv_cols}")

        # Load experimental data
        from ..data import load_experiment_data

        exp_data = load_experiment_data(self.config.exp_csv)
        power_points = exp_data[:, 0]
        exp_observations = exp_data[:, 1:4]
        output_weights = np.asarray(self.config.output_weights, dtype=float)

        # Detect power point columns
        power_cols_map = {}
        for power in power_points:
            w_col = f"width_um_{int(power)}W"
            d_col = f"depth_um_{int(power)}W"
            a_col = f"area_um2_{int(power)}W"
            if all(c in df.columns for c in [w_col, d_col, a_col]):
                power_cols_map[power] = (w_col, d_col, a_col)

        has_detailed_data = len(power_cols_map) == len(power_points)

        for _, row in df.iterrows():
            params = []
            for param_name in self.runner.active_names:
                csv_col = self._resolve_param_csv_column(df, param_name)
                params.append(float(row[csv_col]))

            predictions = None
            nrmse = np.nan

            # Prefer the objective already stored in the CSV (during warm start, historical objective values are not overwritten by recomputation)
            if "objective" in df.columns and not pd.isna(row["objective"]):
                nrmse = float(row["objective"])

            # Attempt to parse detailed data
            if has_detailed_data:
                try:
                    current_preds = []
                    any_nan = False
                    for power in power_points:
                        w_col, d_col, a_col = power_cols_map[power]
                        vals = [row[w_col], row[d_col], row[a_col]]
                        if any(pd.isna(v) for v in vals):
                            any_nan = True
                            break
                        current_preds.append(vals)

                    if not any_nan:
                        predictions = np.array(current_preds)
                        if np.isnan(nrmse):
                            nrmse, _ = compute_bayes_nrmse_percent(
                                predictions,
                                exp_observations,
                                weights=output_weights,
                            )
                except Exception:
                    predictions = None

            # If objective is missing and cannot be recomputed, fall back to legacy NRMSE column name
            if np.isnan(nrmse):
                if "NRMSE" in df.columns and not pd.isna(row["NRMSE"]):
                    nrmse = float(row["NRMSE"])

            all_x.append(params)
            all_y.append(nrmse)
            all_predictions.append(predictions)

        print(f"  Loaded {len(all_x)} records from {csv_path}")
        return all_x, all_y, all_predictions
