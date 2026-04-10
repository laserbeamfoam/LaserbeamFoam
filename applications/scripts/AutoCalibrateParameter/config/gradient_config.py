"""
梯度优化配置

用于 scipy.optimize.least_squares 点估计优化的参数配置
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from .base_config import BaseConfig


@dataclass
class GradientConfig(BaseConfig):
    """
    梯度优化配置 (scipy.optimize.least_squares)

    继承 BaseConfig，添加梯度优化器参数

    Attributes
    ----------
    max_nfev : int
        最大函数评估次数
    ftol : float
        目标函数相对收敛容差
    xtol : float
        参数相对收敛容差
    diff_step : float
        有限差分步长
    k_lower : float
        参数缩放因子下界
    k_upper : float
        参数缩放因子上界
    n_directions : int
        多起点优化的方向数
    init_strategy : str
        初始化策略
    direction_scale : float
        方向缩放因子
    max_workers : int
        并行优化最大进程数
    log_every_eval : int
        每隔多少次评估记录一次日志
    monitor_interval_sec : float
        监控线程轮询间隔 (秒)
    name : str
        优化器名称
    """

    # === 优化参数 ===
    max_nfev: int = 300
    ftol: float = 1e-4
    xtol: float = 1e-4
    diff_step: float = 0.05

    # === 参数缩放因子边界 ===
    k_lower: float = 0.5
    k_upper: float = 1.5

    # === 多起点优化 ===
    n_directions: int = 8
    init_strategy: Literal["legacy_single", "legacy_multi", "directional"] = "legacy_multi"
    direction_scale: float = 0.5

    # === 并行 ===
    max_workers: int = 8

    # === 监控 ===
    log_every_eval: int = 1
    monitor_interval_sec: float = 0.5

    # === 输出 ===
    name: str = "meltpool_gradient"

    def to_dict(self) -> dict:
        """转换为字典"""
        base_dict = super().to_dict()
        base_dict.update({
            "max_nfev": self.max_nfev,
            "ftol": self.ftol,
            "xtol": self.xtol,
            "diff_step": self.diff_step,
            "k_lower": self.k_lower,
            "k_upper": self.k_upper,
            "n_directions": self.n_directions,
            "init_strategy": self.init_strategy,
            "direction_scale": self.direction_scale,
            "max_workers": self.max_workers,
            "log_every_eval": self.log_every_eval,
            "monitor_interval_sec": self.monitor_interval_sec,
            "name": self.name,
        })
        return base_dict
