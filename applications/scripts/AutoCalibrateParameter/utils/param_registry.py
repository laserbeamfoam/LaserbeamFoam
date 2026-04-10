"""
参数注册中心

集中管理参数顺序、默认边界、CSV列名映射与显示标签。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple


@dataclass(frozen=True)
class ParameterSpec:
    """单个参数的元数据定义。"""

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
    """返回按优化顺序排列的参数名称。"""
    return [spec.name for spec in PARAMETER_SPECS]


def get_param_labels() -> Dict[str, str]:
    """返回参数显示标签映射。"""
    return {spec.name: spec.label for spec in PARAMETER_SPECS}


def get_param_latex_labels() -> List[str]:
    """返回 LaTeX 标签顺序列表。"""
    return [spec.latex_label for spec in PARAMETER_SPECS]


def get_default_bounds_map() -> Dict[str, Tuple[float, float]]:
    """返回参数默认边界映射。"""
    return {spec.name: spec.default_bounds for spec in PARAMETER_SPECS}


def get_param_csv_column(name: str) -> str:
    """返回参数在历史 CSV 中对应的列名。"""
    return _PARAM_SPEC_MAP[name].csv_column


def canonical_param_name(name: str) -> str:
    """返回参数规范名（兼容旧别名）。"""
    return _PARAM_SPEC_MAP[name].name


def get_param_csv_candidates(name: str) -> List[str]:
    """返回参数在历史 CSV 中可接受的列名（新列名 + 兼容旧列名）。"""
    spec = _PARAM_SPEC_MAP[name]
    return [spec.csv_column, *spec.legacy_csv_columns]


def get_param_csv_column_map() -> Dict[str, str]:
    """返回参数名 -> CSV列名映射。"""
    return {spec.name: spec.csv_column for spec in PARAMETER_SPECS}


def get_param_bounds_from_config(config) -> Dict[str, Tuple[float, float]]:
    """从配置对象中提取参数边界；缺失时使用注册中心默认值。"""
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
