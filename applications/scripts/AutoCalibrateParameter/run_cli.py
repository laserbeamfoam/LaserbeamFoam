#!/usr/bin/env python3
"""
Melt pool parameter calibration main entry point

Unified command line interface supporting two optimization methods:
- bayes: Bayesian optimization (scikit-optimize)
- gradient: gradient optimization (scipy.optimize.least_squares)

Usage examples:
    # Bayesian optimization
    python main.py --method bayes --n-batches 50 --batch-size 5

    # Gradient optimization
    python main.py --method gradient --n-directions 8 --max-nfev 300
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Optional

import numpy as np


def parse_args() -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Melt pool parameter calibration tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # General parameters
    parser.add_argument(
        "--method",
        type=str,
        choices=["bayes", "gradient"],
        default="bayes",
        help="Optimization method: bayes (Bayesian optimization), gradient (gradient optimization) (default: bayes)",
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config.yaml",
        help="Configuration file path (default: config.yaml)",
    )
    parser.add_argument(
        "--exp-csv",
        type=str,
        default=None,
        help="Experimental data CSV file path (overrides config file)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Output directory (overrides config file)",
    )
    parser.add_argument(
        "--n-proc",
        type=int,
        default=None,
        help="Number of parallel OpenFOAM processes (overrides config file)",
    )
    parser.add_argument(
        "--hpc",
        action="store_true",
        help="Enable HPC mode",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Verbose output",
    )

    # Bayesian optimization specific parameters
    bayes_group = parser.add_argument_group("Bayesian optimization parameters")
    bayes_group.add_argument(
        "--n-initial-points",
        type=int,
        default=None,
        help="Number of initial sampling points (default: 15)",
    )
    bayes_group.add_argument(
        "--n-batches",
        type=int,
        default=None,
        help="Number of optimization batches (default: 50)",
    )
    bayes_group.add_argument(
        "--batch-size",
        type=int,
        default=None,
        help="Number of evaluations per batch (default: 5)",
    )
    bayes_group.add_argument(
        "--acq-func",
        type=str,
        choices=["EI", "PI", "LCB", "gp_hedge"],
        default=None,
        help="Acquisition function (default: EI)",
    )
    bayes_group.add_argument(
        "--warmstart-csv",
        type=str,
        default=None,
        help="Historical data CSV file path for warm start (e.g. bayes_history.csv)",
    )

    # Gradient optimization specific parameters
    gradient_group = parser.add_argument_group("Gradient optimization parameters")
    gradient_group.add_argument(
        "--n-directions",
        type=int,
        default=None,
        help="Number of multi-start directions (default: 8)",
    )
    gradient_group.add_argument(
        "--max-nfev",
        type=int,
        default=None,
        help="Maximum number of function evaluations (default: 300)",
    )
    gradient_group.add_argument(
        "--ftol",
        type=float,
        default=None,
        help="Objective function convergence tolerance (default: 1e-4)",
    )
    gradient_group.add_argument(
        "--xtol",
        type=float,
        default=None,
        help="Parameter convergence tolerance (default: 1e-4)",
    )
    gradient_group.add_argument(
        "--init-strategy",
        type=str,
        choices=["legacy_single", "legacy_multi", "directional"],
        default=None,
        help="Initialization strategy (default: legacy_multi)",
    )

    return parser.parse_args()


def create_config(args: argparse.Namespace):
    """Create a configuration object from command line arguments"""
    from AutoCalibrateParameter.config import BayesConfig, GradientConfig

    # Base paths
    script_dir = Path(__file__).parent
    base_dir = script_dir

    # Check if config exists in CWD or absolute path
    if Path(args.config).exists():
        config_path = Path(args.config).resolve()
    else:
        # Fallback to script directory (but we prefer CWD for runs)
        config_path = script_dir / args.config

    # Default runs root to CWD/runs if not specified in config
    default_runs_root = Path.cwd() / "runs"

    # Select config class based on method
    if args.method == "bayes":
        if config_path.exists():
            config = BayesConfig.from_yaml(config_path)
        else:
            config = BayesConfig(
                case_dir=base_dir / "openfoam_case",
                postproc_script=base_dir / "postprocess.py",
                exp_csv=base_dir / "experimental_data.csv",
                runs_root=default_runs_root,
            )

        # Override with command line arguments
        if args.n_initial_points is not None:
            config.n_initial_points = args.n_initial_points
        if args.n_batches is not None:
            config.n_batches = args.n_batches
        if args.batch_size is not None:
            config.batch_size = args.batch_size
        if args.acq_func is not None:
            config.acq_func = args.acq_func
        if args.warmstart_csv is not None:
            config.warmstart_csv = args.warmstart_csv

    elif args.method == "gradient":
        if config_path.exists():
            config = GradientConfig.from_yaml(config_path)
        else:
            config = GradientConfig(
                case_dir=base_dir / "openfoam_case",
                postproc_script=base_dir / "postprocess.py",
                exp_csv=base_dir / "experimental_data.csv",
                runs_root=base_dir / "runs",
            )

        # Override with command line arguments
        if args.n_directions is not None:
            config.n_directions = args.n_directions
        if args.max_nfev is not None:
            config.max_nfev = args.max_nfev
        if args.ftol is not None:
            config.ftol = args.ftol
        if args.xtol is not None:
            config.xtol = args.xtol
        if args.init_strategy is not None:
            config.init_strategy = args.init_strategy

    else:
        raise ValueError(f"Unknown method: {args.method}")

    # General parameter overrides
    if args.exp_csv is not None:
        config.exp_csv = Path(args.exp_csv)
    if args.output_dir is not None:
        config.runs_root = Path(args.output_dir)
    if args.n_proc is not None:
        config.n_proc = args.n_proc
    if args.hpc:
        config.hpc_mode = True

    return config


def run_bayes(config, exp_data: np.ndarray, verbose: bool = False):
    """Run Bayesian optimization"""
    from AutoCalibrateParameter.optimizers import BayesOptimizer
    from AutoCalibrateParameter.simulation import SimulationRunner

    print("\n" + "=" * 60)
    print("Bayesian optimization (scikit-optimize)")
    print("=" * 60)
    print(f"Initial points: {config.n_initial_points}")
    print(f"Number of batches: {config.n_batches}")
    print(f"Batch size: {config.batch_size}")
    print(f"Acquisition function: {config.acq_func}")
    print()

    runner = SimulationRunner(config)
    optimizer = BayesOptimizer(config, runner)
    result = optimizer.optimize(exp_data)

    return result


def run_gradient(config, exp_data: np.ndarray, verbose: bool = False):
    """Run gradient optimization"""
    from AutoCalibrateParameter.optimizers import GradientOptimizer
    from AutoCalibrateParameter.simulation import SimulationRunner

    print("\n" + "=" * 60)
    print("Gradient optimization (scipy.optimize.least_squares)")
    print("=" * 60)
    print(f"Number of directions: {config.n_directions}")
    print(f"Max evaluations: {config.max_nfev}")
    print(f"Initialization strategy: {config.init_strategy}")
    print()

    runner = SimulationRunner(config)
    optimizer = GradientOptimizer(config, runner)
    result = optimizer.optimize(exp_data)

    return result


def print_result(result, param_names: list = None):
    """Print optimization results"""
    if param_names is None:
        from AutoCalibrateParameter.utils import get_param_names

        param_names = get_param_names()

    print("\n" + "=" * 60)
    print("Optimization Results")
    print("=" * 60)
    print(f"Method: {result.method}")
    print(f"Best objective value: {result.best_cost:.6e}")
    print(f"Total evaluations: {result.n_evaluations}")
    print(f"Message: {result.message}")
    print()
    print("Optimal parameters:")
    for name, val in zip(param_names, result.best_params):
        print(f"  {name}: {val:.6g}")


def main():
    """Main function"""
    args = parse_args()

    print("=" * 60)
    print("Melt pool parameter calibration tool")
    print("=" * 60)
    print(f"Method: {args.method}")
    print(f"Config: {args.config}")
    print()

    # Create configuration
    try:
        config = create_config(args)
    except Exception as e:
        print(f"Error: cannot create configuration - {e}")
        sys.exit(1)

    # Load experimental data
    from AutoCalibrateParameter.data import load_experiment_data

    try:
        exp_data = load_experiment_data(config.exp_csv)
        print(f"Loading experimental data: {config.exp_csv}")
        print(f"  Number of data points: {len(exp_data)}")
        print(f"  Power range: {exp_data[:, 0].min():.0f} - {exp_data[:, 0].max():.0f} W")
    except Exception as e:
        print(f"Error: cannot load experimental data - {e}")
        sys.exit(1)

    # Create output directory
    # Use runs_root directly; no longer auto-append "{method}_optimization" subdirectory
    output_dir = config.runs_root

    # Compatibility with legacy configs: if runs_root already looks like ".../bayes_optimization",
    # fall back to the parent runs directory
    legacy_method_dir = f"{args.method}_optimization"
    if output_dir.name == legacy_method_dir:
        output_dir = output_dir.parent

    output_dir.mkdir(parents=True, exist_ok=True)

    # Important: update config.runs_root to ensure CSV is saved to the correct location
    config.runs_root = output_dir
    print(f"Output directory: {output_dir}")

    # Save configuration
    config_save_path = output_dir / "config.yaml"
    config.save_yaml(config_save_path)
    print(f"Config saved: {config_save_path}")

    # Run optimization
    start_time = time.time()

    try:
        if args.method == "bayes":
            result = run_bayes(config, exp_data, args.verbose)
        elif args.method == "gradient":
            result = run_gradient(config, exp_data, args.verbose)
        else:
            raise ValueError(f"Unknown method: {args.method}")
    except KeyboardInterrupt:
        print("\nUser interrupt")
        sys.exit(1)
    except Exception as e:
        print(f"\nError: optimization failed - {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    elapsed = time.time() - start_time
    print(f"\nTotal elapsed time: {elapsed / 60:.1f} minutes")

    # Print results
    print_result(result)

    # Save results
    result_path = output_dir / "result.json"
    result.save(result_path)
    print(f"\nResults saved: {result_path}")

    # Generate comparison plots
    try:
        from AutoCalibrateParameter.postprocess import (
            save_comparison_plot,
            plot_predictions_vs_experiments,
        )
        from AutoCalibrateParameter.simulation import SimulationRunner

        runner = SimulationRunner(config)
        power_points = exp_data[:, 0]
        observations = exp_data[:, 1:4]
        predictions = runner.run(result.best_params, power_points)

        plots_dir = output_dir / "plots"
        plots_dir.mkdir(exist_ok=True)

        save_comparison_plot(
            plots_dir,
            power_points,
            observations,
            predictions,
            nrmse=result.best_cost,
            title=f"{result.method} optimization results",
        )

        plot_predictions_vs_experiments(
            plots_dir,
            observations,
            predictions,
            title=f"{result.method} predictions vs experiments",
        )

        print(f"Plots saved: {plots_dir}")

    except Exception as e:
        print(f"Warning: cannot generate plots - {e}")

    print("\nDone!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
