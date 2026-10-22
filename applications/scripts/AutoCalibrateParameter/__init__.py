"""
AutoCalibrateParameter - LaserbeamFoam melt pool simulation parameter calibration library

Supports two optimization methods:
- Bayesian: Bayesian optimization (skopt point estimation)
- Gradient: gradient optimization (scipy point estimation)
"""

__version__ = "1.0.0"
__author__ = "CGH"

# Lazy imports — heavy dependencies (skopt, sklearn, scipy) are only loaded
# when the corresponding name is actually accessed, keeping CLI startup fast.

__all__ = [
    # Config
    "BaseConfig",
    "BayesConfig",
    "GradientConfig",
    # Optimizers
    "BaseOptimizer",
    "OptimizationResult",
    "BayesOptimizer",
    "GradientOptimizer",
    # Simulation
    "SimulationRunner",
    # Models
    "OpenFOAMCaseManager",
    # Data
    "load_experiment_data",
]

_LAZY_IMPORTS = {
    "BaseConfig": ".config",
    "BayesConfig": ".config",
    "GradientConfig": ".config",
    "BaseOptimizer": ".optimizers",
    "OptimizationResult": ".optimizers",
    "BayesOptimizer": ".optimizers",
    "GradientOptimizer": ".optimizers",
    "SimulationRunner": ".simulation",
    "OpenFOAMCaseManager": ".models",
    "load_experiment_data": ".data",
}


def __getattr__(name):
    if name in _LAZY_IMPORTS:
        import importlib
        module = importlib.import_module(_LAZY_IMPORTS[name], __name__)
        value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
