"""
贝叶斯优化配置

用于 scikit-optimize 点估计优化的参数配置
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from .base_config import BaseConfig


@dataclass
class BayesConfig(BaseConfig):
    """
    贝叶斯优化配置 (scikit-optimize)

    继承 BaseConfig，添加 skopt 优化器参数

    Attributes
    ----------
    n_initial_points : int
        初始随机采样点数
    n_batches : int
        贝叶斯优化批次数
    batch_size : int
        每批次并行评估的候选点数
    acq_func : str
        采集函数 ("EI", "PI", "LCB", "gp_hedge")
    xi : float
        EI/PI 的探索参数
    kappa : float
        LCB 的探索参数
    noise : float
        GP 噪声参数
    k_lower : float
        参数缩放因子下界
    k_upper : float
        参数缩放因子上界
    max_workers : int
        并行评估最大进程数
    resume : bool
        是否从上次中断处恢复
    name : str
        优化器名称
    """

    # === 采样配置 ===
    n_initial_points: int = 15
    n_batches: int = 1
    batch_size: int = 1

    # === 采集函数 ===
    acq_func: Literal["EI", "PI", "LCB", "gp_hedge"] = "EI"
    xi: float = 0.01
    kappa: float = 1.96
    noise: float = 1e-10

    # === 参数缩放因子边界 ===
    k_lower: float = 0.5
    k_upper: float = 1.5

    # === 并行与恢复 ===
    max_workers: int = 1
    resume: bool = True
    warmstart_csv: str = "bayes_history.csv"  # 历史数据 CSV 文件路径，用于热启动

    # === 输出 ===
    name: str = "meltpool_bayes"

    def to_dict(self) -> dict:
        """转换为字典"""
        base_dict = super().to_dict()
        base_dict.update({
            "n_initial_points": self.n_initial_points,
            "n_batches": self.n_batches,
            "batch_size": self.batch_size,
            "acq_func": self.acq_func,
            "xi": self.xi,
            "kappa": self.kappa,
            "noise": self.noise,
            "k_lower": self.k_lower,
            "k_upper": self.k_upper,
            "max_workers": self.max_workers,
            "resume": self.resume,
            "warmstart_csv": self.warmstart_csv,
            "name": self.name,
        })
        return base_dict
