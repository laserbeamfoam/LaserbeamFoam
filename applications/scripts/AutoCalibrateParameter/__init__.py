"""
AutoCalibrateParameter - LaserbeamFoam 熔池仿真参数校准库

支持两种优化方法:
- Bayesian: 贝叶斯优化 (skopt 点估计)
- Gradient: 梯度优化 (scipy 点估计)
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
