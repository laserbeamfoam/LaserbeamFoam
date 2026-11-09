"""
Base configuration class

Configuration parameters shared by all optimization methods
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
import re
from typing import Any, Dict, List, Optional, Tuple, Union, get_args, get_origin, get_type_hints

import yaml

from ..utils.param_registry import canonical_param_name

# Package root directory (applications/scripts/AutoCalibrateParameter)
_PACKAGE_ROOT = Path(__file__).resolve().parent.parent
# Repository root directory (LaserbeamFoam)
REPO_ROOT = _PACKAGE_ROOT.parent.parent.parent
# Project case root directory (Project/AutoCalibrateParameter)
PROJECT_ROOT = REPO_ROOT / "Project" / "AutoCalibrateParameter"

_INT_LITERAL_RE = re.compile(r"^[+-]?\d+$")
_NUMERIC_LITERAL_RE = re.compile(
    r"^[+-]?(?:(?:\d+\.\d*)|(?:\d+)|(?:\.\d+))(?:[eE][+-]?\d+)?$"
)


def _parse_numeric_literal(value: str):
    """Try parsing plain/scientific numeric strings; return original on failure."""
    if not isinstance(value, str):
        return value

    text = value.strip()
    if not text or not _NUMERIC_LITERAL_RE.match(text):
        return value

    if _INT_LITERAL_RE.match(text):
        try:
            return int(text)
        except ValueError:
            return value

    try:
        return float(text)
    except ValueError:
        return value


def _coerce_typed_value(value: Any, expected_type: Any):
    """Coerce YAML-loaded values to declared dataclass field types when safe."""
    origin = get_origin(expected_type)
    args = get_args(expected_type)

    if expected_type is float:
        parsed = _parse_numeric_literal(value)
        if isinstance(parsed, bool):
            return value
        if isinstance(parsed, (int, float)):
            return float(parsed)
        return value

    if expected_type is int:
        parsed = _parse_numeric_literal(value)
        if isinstance(parsed, bool):
            return value
        if isinstance(parsed, int):
            return parsed
        if isinstance(parsed, float) and parsed.is_integer():
            return int(parsed)
        return value

    if origin in (list, List):
        if not isinstance(value, list):
            return value
        if not args:
            return value
        return [_coerce_typed_value(item, args[0]) for item in value]

    if origin in (tuple, Tuple):
        if isinstance(value, list):
            value = tuple(value)
        if not isinstance(value, tuple):
            return value
        if not args:
            return value

        if len(args) == 2 and args[1] is Ellipsis:
            return tuple(_coerce_typed_value(item, args[0]) for item in value)

        coerced = []
        for index, item in enumerate(value):
            item_type = args[index] if index < len(args) else args[-1]
            coerced.append(_coerce_typed_value(item, item_type))
        return tuple(coerced)

    if origin in (dict, Dict):
        if not isinstance(value, dict):
            return value
        if len(args) != 2:
            return value
        key_type, value_type = args
        return {
            _coerce_typed_value(key, key_type): _coerce_typed_value(val, value_type)
            for key, val in value.items()
        }

    if origin is Union:
        if value is None:
            return value
        non_none_options = [option for option in args if option is not type(None)]
        for option in non_none_options:
            coerced = _coerce_typed_value(value, option)
            if coerced is not value:
                return coerced
        return value

    return value


def _default_case_dir() -> Path:
    return PROJECT_ROOT


def _default_postproc_script() -> Path:
    return REPO_ROOT / "applications" / "scripts" / "postProcessing" / "characterise_meltpool.py"


def _default_exp_csv() -> Path:
    return PROJECT_ROOT / "experimental_data.csv"


def _default_runs_root() -> Path:
    return PROJECT_ROOT / "runs"


def _registry_param_names_fallback() -> List[str]:
    """Retrieve parameter names from the parameter registry; fall back to default order on failure."""
    fallback = ["sigma", "marangoni", "substrate_temp", "absorptivity", "recoilCoeff", "radius_flavour"]
    try:
        from ..utils.param_registry import PARAM_NAMES

        return list(PARAM_NAMES)
    except Exception:
        return fallback


def _registry_bounds_from_config(config: "BaseConfig") -> Dict[str, Tuple[float, float]]:
    """Extract parameter bounds from the config object via the registry; fall back to explicit fields on failure."""
    try:
        from ..utils.param_registry import get_param_bounds_from_config

        return get_param_bounds_from_config(config)
    except Exception:
        return {
            "sigma": tuple(config.sigma_bounds),
            "marangoni": tuple(config.marangoni_bounds),
            "substrate_temp": tuple(config.substrate_temp_bounds),
            "absorptivity": tuple(config.absorptivity_bounds),
            "recoilCoeff": tuple(
                getattr(config, "recoilCoeff_bounds", getattr(config, "damper_bounds"))
            ),
            "radius_flavour": tuple(config.radius_flavour_bounds),
        }


@dataclass
class BaseConfig:
    """
    Base configuration shared by all optimization methods

    Attributes
    ----------
    case_dir : Path
        OpenFOAM case directory
    postproc_script : Path
        Post-processing script path
    exp_csv : Path
        Experimental data CSV file path
    runs_root : Path
        Root directory for run results
    foam_bashrc : str, optional
        OpenFOAM bashrc path
    foam_runner : str, optional
        OpenFOAM runner command (e.g. "of2506")
    n_proc : int
        Number of MPI parallel cores
    hpc_mode : bool
        Whether to use HPC mode
    mpirun_flags : tuple
        Additional mpirun arguments
    postproc_python : str
        Post-processing Python interpreter
    postproc_runner : str, optional
        Post-processing runner
    pvpython : str, optional
        pvpython path
    sigma_bounds : tuple
        sigma parameter bounds
    marangoni_bounds : tuple
        Marangoni constant bounds
    substrate_temp_bounds : tuple
        Substrate temperature bounds (K)
    absorptivity_bounds : tuple
        Absorptivity bounds
    recoilCoeff_bounds : tuple
        Recoil pressure coefficient bounds
    radius_flavour_bounds : tuple
        Spot radial distribution factor Radius_Flavour bounds
    power_range : tuple
        Laser power range (W)
    seed : int, optional
        Random seed (for reproducibility)
    laser_diameter : float
        Laser diameter (m), used for post-processing
    cell_size : float
        Mesh cell size (m), used for post-processing
    x_domain : tuple
        Computational domain X-direction range (m)
    y_begin_track : float
        Track analysis start Y coordinate (m)
    y_end_track : float
        Track analysis end Y coordinate (m)
    plot_geometry : bool
        Whether to generate post-processing plots
    archive_mode : str
        Result archiving mode: "latest" saves only the last frame; "all" saves all time-step results
    """

    # === Path configuration ===
    case_dir: Path = field(default_factory=_default_case_dir)
    postproc_script: Path = field(default_factory=_default_postproc_script)
    exp_csv: Path = field(default_factory=_default_exp_csv)
    runs_root: Path = field(default_factory=_default_runs_root)

    # === OpenFOAM configuration ===
    foam_bashrc: Optional[str] = None
    foam_runner: Optional[str] = None
    n_proc: int = 12
    hpc_mode: bool = False
    mpirun_flags: Tuple[str, ...] = ("--oversubscribe",)
    postproc_python: str = "python"
    postproc_runner: Optional[str] = None
    pvpython: Optional[str] = None

    # === Parameter bounds ===
    sigma_bounds: Tuple[float, float] = (1.0, 2.0)
    marangoni_bounds: Tuple[float, float] = (-8e-4, -4e-6)
    substrate_temp_bounds: Tuple[float, float] = (300.0, 800.0)
    absorptivity_bounds: Tuple[float, float] = (0.5, 3.0)
    recoilCoeff_bounds: Tuple[float, float] = (0.5, 2.0)
    damper_bounds: Tuple[float, float] = (0.5, 2.0)  # compatibility with legacy config
    radius_flavour_bounds: Tuple[float, float] = (1.0, 3.0)
    laser_radius_bounds: Tuple[float, float] = (30e-6, 70e-6)

    # === Power range ===
    power_range: Tuple[float, float] = (140.0, 260.0)

    # === Output weight configuration ===
    output_weights: List[float] = field(default_factory=lambda: [1.0, 1.0, 0.0])

    # === Selective parameter optimization configuration ===
    active_params: Optional[List[str]] = None  # None means optimize all parameters
    fixed_values: Dict[str, float] = field(default_factory=dict)

    # === Post-processing geometry parameters ===
    laser_diameter: float = 70e-6
    cell_size: float = 1e-6
    x_domain: Tuple[float, float] = (0.0, 0.0003)
    y_begin_track: float = 150e-6
    y_end_track: float = 300e-6
    plot_geometry: bool = True

    # === Result archiving configuration ===
    archive_mode: str = "latest"

    # === Reproducibility ===
    seed: Optional[int] = 42

    def __post_init__(self):
        """Ensure correct path types and validate configuration"""
        self.case_dir = Path(self.case_dir)
        self.postproc_script = Path(self.postproc_script)
        self.exp_csv = Path(self.exp_csv)
        self.runs_root = Path(self.runs_root)

        # Validate output_weights length
        if len(self.output_weights) != 3:
            raise ValueError("output_weights must have 3 elements [width, depth, area]")

        # Validate archive mode
        archive_mode = "latest" if self.archive_mode is None else str(self.archive_mode).strip().lower()
        if archive_mode not in {"latest", "all"}:
            raise ValueError("archive_mode must be 'latest' or 'all'")
        self.archive_mode = archive_mode

        all_param_names = _registry_param_names_fallback()

        # Validate active_params
        if self.active_params is not None:
            canonical_active = [canonical_param_name(name) for name in self.active_params]
            invalid = set(canonical_active) - set(all_param_names)
            if invalid:
                raise ValueError(f"Invalid active_params: {invalid}")
            self.active_params = canonical_active

        # Validate fixed_values
        if self.fixed_values:
            for name in self.fixed_values:
                canonical_name = canonical_param_name(name)
                if canonical_name not in all_param_names:
                    raise ValueError(f"Unknown parameter in fixed_values: {name}")
            self.fixed_values = {
                canonical_param_name(name): value
                for name, value in self.fixed_values.items()
            }

    def get_param_bounds(self) -> List[Tuple[float, float]]:
        """
        Return the list of parameter bounds

        Returns
        -------
        list of tuple
            [(sigma_min, sigma_max), (marangoni_min, marangoni_max), ...]
        """
        bounds_map = _registry_bounds_from_config(self)
        return [bounds_map[name] for name in _registry_param_names_fallback()]

    def get_param_names(self) -> List[str]:
        """Return the list of parameter names"""
        return _registry_param_names_fallback()

    def build_runs_path(self, method: str) -> Path:
        """
        Build a timestamped run directory

        Parameters
        ----------
        method : str
            Optimization method name (bayes, gradient)

        Returns
        -------
        Path
            Path of the form runs/runs_bayes_20250119_120000
        """
        safe_method = method.strip().lower().replace(" ", "_")
        if not safe_method:
            safe_method = "unknown"
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return self.runs_root / f"runs_{safe_method}_{timestamp}"

    @classmethod
    def from_yaml(cls, yaml_path: Path) -> "BaseConfig":
        """
        Load configuration from a YAML file

        Parameters
        ----------
        yaml_path : Path
            Path to the YAML configuration file

        Returns
        -------
        BaseConfig
            Configuration instance
        """
        yaml_path = Path(yaml_path)
        if not yaml_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {yaml_path}")

        with open(yaml_path, "r") as f:
            data = yaml.safe_load(f) or {}

        if not isinstance(data, dict):
            raise ValueError(f"Configuration file format error, expected a dictionary: {yaml_path}")

        # Strict validation: disallow unknown fields
        valid_fields = {f.name for f in cls.__dataclass_fields__.values()}
        unknown_fields = sorted(set(data.keys()) - valid_fields)
        if unknown_fields:
            raise ValueError(f"Configuration file contains unknown fields: {unknown_fields}")

        # Type-safe numeric conversion based on dataclass annotations (e.g. 70e-6 -> float)
        type_hints = get_type_hints(cls)
        for field_name, expected_type in type_hints.items():
            if field_name in data:
                data[field_name] = _coerce_typed_value(data[field_name], expected_type)

        # Handle paths (relative paths are resolved relative to the config file directory)
        config_dir = yaml_path.parent
        path_fields = {"case_dir", "postproc_script", "exp_csv", "runs_root"}
        for field_name in path_fields:
            if field_name in data:
                field_path = Path(data[field_name])
                if not field_path.is_absolute():
                    data[field_name] = config_dir / field_path

        # Convert mpirun_flags to tuple
        if "mpirun_flags" in data and isinstance(data["mpirun_flags"], list):
            data["mpirun_flags"] = tuple(data["mpirun_flags"])

        # Convert x_domain to tuple
        if "x_domain" in data and isinstance(data["x_domain"], list):
            data["x_domain"] = tuple(data["x_domain"])

        return cls(**data)

    def to_dict(self) -> dict:
        """Convert to dictionary (for serialization)"""
        return {
            "case_dir": str(self.case_dir),
            "postproc_script": str(self.postproc_script),
            "exp_csv": str(self.exp_csv),
            "runs_root": str(self.runs_root),
            "foam_bashrc": self.foam_bashrc,
            "foam_runner": self.foam_runner,
            "n_proc": self.n_proc,
            "hpc_mode": self.hpc_mode,
            "mpirun_flags": list(self.mpirun_flags),
            "postproc_python": self.postproc_python,
            "postproc_runner": self.postproc_runner,
            "pvpython": self.pvpython,
            "sigma_bounds": list(self.sigma_bounds),
            "marangoni_bounds": list(self.marangoni_bounds),
            "substrate_temp_bounds": list(self.substrate_temp_bounds),
            "absorptivity_bounds": list(self.absorptivity_bounds),
            "recoilCoeff_bounds": list(self.recoilCoeff_bounds),
            "radius_flavour_bounds": list(self.radius_flavour_bounds),
            "laser_radius_bounds": list(self.laser_radius_bounds),
            "power_range": list(self.power_range),
            "output_weights": self.output_weights,
            "active_params": self.active_params,
            "fixed_values": self.fixed_values,
            "seed": self.seed,
            "laser_diameter": self.laser_diameter,
            "cell_size": self.cell_size,
            "x_domain": list(self.x_domain),
            "y_begin_track": self.y_begin_track,
            "y_end_track": self.y_end_track,
            "plot_geometry": self.plot_geometry,
            "archive_mode": self.archive_mode,
        }

    def save_yaml(self, yaml_path: Path) -> None:
        """
        Save configuration to a YAML file

        Parameters
        ----------
        yaml_path : Path
            Path to save the YAML file
        """
        yaml_path = Path(yaml_path)
        yaml_path.parent.mkdir(parents=True, exist_ok=True)

        with open(yaml_path, "w") as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False, allow_unicode=True)
