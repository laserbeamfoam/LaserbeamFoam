"""
Parameter registry

Centrally manages parameter order, default bounds, CSV column name mapping, and display labels.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple


@dataclass(frozen=True)
class ParameterSpec:
    """Metadata definition for a single parameter."""

    name: str
    label: str
    latex_label: str
    default_bounds: Tuple[float, float]
    config_bounds_attr: str
    csv_column: str
    legacy_names: Tuple[str, ...] = field(default_factory=tuple)
    legacy_bounds_attrs: Tuple[str, ...] = field(default_factory=tuple)
    legacy_csv_columns: Tuple[str, ...] = field(default_factory=tuple)


PARAMETER_SPECS: Tuple[ParameterSpec, ...] = (
    ParameterSpec(
        name="sigma",
        label=r"$\sigma$ (N/m)",
        latex_label=r"$\sigma$",
        default_bounds=(1.0, 2.0),
        config_bounds_attr="sigma_bounds",
        csv_column="sigma",
    ),
    ParameterSpec(
        name="marangoni",
        label=r"$\gamma$ (N/m·K)",
        latex_label=r"$\gamma$",
        default_bounds=(-8e-4, -4e-6),
        config_bounds_attr="marangoni_bounds",
        csv_column="Marangoni_Constant",
    ),
    ParameterSpec(
        name="substrate_temp",
        label=r"$T_s$ (K)",
        latex_label=r"$T_s$",
        default_bounds=(300.0, 800.0),
        config_bounds_attr="substrate_temp_bounds",
        csv_column="substrate_temp",
    ),
    ParameterSpec(
        name="absorptivity",
        label=r"$\eta$",
        latex_label=r"$\eta$",
        default_bounds=(0.5, 3.0),
        config_bounds_attr="absorptivity_bounds",
        csv_column="absorptivity",
    ),
    ParameterSpec(
        name="recoilCoeff",
        label=r"$\alpha_{recoil}$",
        latex_label=r"$\alpha_{recoil}$",
        default_bounds=(0.5, 2.0),
        config_bounds_attr="recoilCoeff_bounds",
        csv_column="recoilCoeff",
        legacy_names=("damper",),
        legacy_bounds_attrs=("damper_bounds",),
        legacy_csv_columns=("damper",),
    ),
    ParameterSpec(
        name="radius_flavour",
        label="Radius_Flavour",
        latex_label=r"$F_r$",
        default_bounds=(1.0, 3.0),
        config_bounds_attr="radius_flavour_bounds",
        csv_column="radius_flavour",
    ),
    ParameterSpec(
        name="laser_radius",
        label="laserRadius (m)",
        latex_label=r"$R_L$",
        default_bounds=(30e-6, 70e-6),
        config_bounds_attr="laser_radius_bounds",
        csv_column="laser_radius",
    ),
)

_PARAM_SPEC_MAP: Dict[str, ParameterSpec] = {}
for _spec in PARAMETER_SPECS:
    _PARAM_SPEC_MAP[_spec.name] = _spec
    for _legacy_name in _spec.legacy_names:
        _PARAM_SPEC_MAP[_legacy_name] = _spec


def get_param_names() -> List[str]:
    """Return parameter names ordered by optimization sequence."""
    return [spec.name for spec in PARAMETER_SPECS]


def get_param_labels() -> Dict[str, str]:
    """Return the parameter display label mapping."""
    return {spec.name: spec.label for spec in PARAMETER_SPECS}


def get_param_latex_labels() -> List[str]:
    """Return the ordered list of LaTeX labels."""
    return [spec.latex_label for spec in PARAMETER_SPECS]


def get_default_bounds_map() -> Dict[str, Tuple[float, float]]:
    """Return the parameter default bounds mapping."""
    return {spec.name: spec.default_bounds for spec in PARAMETER_SPECS}


def get_param_csv_column(name: str) -> str:
    """Return the column name corresponding to the parameter in the history CSV."""
    return _PARAM_SPEC_MAP[name].csv_column


def canonical_param_name(name: str) -> str:
    """Return the canonical parameter name (compatible with legacy aliases)."""
    return _PARAM_SPEC_MAP[name].name


def get_param_csv_candidates(name: str) -> List[str]:
    """Return acceptable column names for the parameter in the history CSV (new name + legacy compatible names)."""
    spec = _PARAM_SPEC_MAP[name]
    return [spec.csv_column, *spec.legacy_csv_columns]


def get_param_csv_column_map() -> Dict[str, str]:
    """Return the parameter name -> CSV column name mapping."""
    return {spec.name: spec.csv_column for spec in PARAMETER_SPECS}


def get_param_bounds_from_config(config) -> Dict[str, Tuple[float, float]]:
    """Extract parameter bounds from the config object; use registry defaults when missing."""
    bounds = {}
    for spec in PARAMETER_SPECS:
        value = None
        for attr in (spec.config_bounds_attr, *spec.legacy_bounds_attrs):
            if hasattr(config, attr):
                candidate = getattr(config, attr)
                if candidate is not None:
                    value = candidate
                    break
        if value is None:
            value = spec.default_bounds
        bounds[spec.name] = tuple(value)
    return bounds


PARAM_NAMES: List[str] = get_param_names()
PARAM_LABELS: Dict[str, str] = get_param_labels()
PARAM_LATEX_LABELS: List[str] = get_param_latex_labels()
DEFAULT_BOUNDS: Dict[str, Tuple[float, float]] = get_default_bounds_map()


__all__ = [
    "ParameterSpec",
    "PARAMETER_SPECS",
    "PARAM_NAMES",
    "PARAM_LABELS",
    "PARAM_LATEX_LABELS",
    "DEFAULT_BOUNDS",
    "get_param_names",
    "get_param_labels",
    "get_param_latex_labels",
    "get_default_bounds_map",
    "get_param_csv_column",
    "canonical_param_name",
    "get_param_csv_candidates",
    "get_param_csv_column_map",
    "get_param_bounds_from_config",
]
