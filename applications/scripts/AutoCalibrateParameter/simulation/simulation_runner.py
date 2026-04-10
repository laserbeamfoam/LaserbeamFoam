"""
仿真运行器

封装 OpenFOAMCaseManager，提供简化的接口用于优化循环
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
    统一仿真运行器

    封装 OpenFOAMCaseManager，提供简化的接口用于优化循环

    Attributes
    ----------
    config : BaseConfig
        配置对象
    case_manager : OpenFOAMCaseManager
        OpenFOAM 案例管理器
    work_dir : Path
        工作目录
    eval_count : int
        评估计数器
    """

    def __init__(
        self,
        config: "BaseConfig",
        work_dir: Optional[Path] = None,
    ):
        """
        初始化仿真运行器

        Parameters
        ----------
        config : BaseConfig
            配置对象
        work_dir : Path, optional
            工作目录，默认为 config.runs_root / "simulations"
        """
        self.config = config
        self.work_dir = work_dir or config.runs_root / "simulations"
        self.work_dir = Path(self.work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)

        # 初始化 OpenFOAM 案例管理器
        self.case_manager = OpenFOAMCaseManager.from_config(config)

        # 参数过滤配置
        self.all_param_names = PARAM_NAMES
        config_bounds = get_param_bounds_from_config(config)

        self.active_names, self.active_bounds, self.fixed_dict = filter_params_for_optimization(
            PARAM_NAMES,
            config_bounds,
            config.active_params,
            config.fixed_values,
        )

        print(f"优化参数: {self.active_names}")
        if self.fixed_dict:
            print(f"固定参数: {self.fixed_dict}")

        # 评估计数器
        self.eval_count = 0

    def get_active_bounds(self) -> List[Tuple[float, float]]:
        """返回优化参数的边界（供优化器使用）"""
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
        运行仿真并返回预测值

        Parameters
        ----------
        params : np.ndarray, shape (n_active,)
            优化参数值（只包含active参数）
        power_points : np.ndarray
            功率点列表 [P1, P2, ...]
        job_id : str | int, optional
            任务 ID，用于归档结果
        on_power_complete : Callable[[int, float, np.ndarray], None], optional
            每个功率点完成后的回调函数
            参数: (power_index, power_value, predictions_so_far)
        start_index : int, optional
            起始功率点索引（用于断点续传），默认为 0

        Returns
        -------
        predictions : np.ndarray, shape (n_powers, 3)
            预测值 [width, depth, area] 单位: μm, μm, μm²
        """
        params = np.asarray(params).flatten()
        power_points = np.asarray(power_points).flatten()

        # 合并active + fixed → 完整参数向量
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
                # 更新参数
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

                # 运行仿真
                self.case_manager.run_simulation()
                self.case_manager.run_postprocess()

                # 读取结果 (m → μm)
                metrics = self.case_manager.read_metrics()
                predictions[i, 0] = metrics["width_mean_m"] * 1e6
                predictions[i, 1] = metrics["depth_mean_m"] * 1e6
                predictions[i, 2] = metrics["area_mean_m2"] * 1e12

                print(
                    f"  [eval {self.eval_count}] P={power:.0f}W: "
                    f"w={predictions[i, 0]:.2f}μm, d={predictions[i, 1]:.2f}μm"
                )

            except Exception as e:
                print(f"  [eval {self.eval_count}] P={power:.0f}W 仿真失败: {e}")
                predictions[i, :] = np.nan

            finally:
                # 归档结果
                if job_id is not None:
                    archive_dir = self.config.runs_root / f"{job_id}" / f"{int(power)}W"
                    archive_dir.mkdir(parents=True, exist_ok=True)
                    archive_mode = getattr(self.config, "archive_mode", "latest")
                    self.case_manager.archive_results(archive_dir, mode=archive_mode)

                # 调用回调函数（如果提供）
                if on_power_complete is not None:
                    try:
                        global_idx = start_index + i
                        on_power_complete(global_idx, power, predictions[: i + 1].copy())
                    except Exception as e:
                        print(f"  [警告] 功率点完成回调失败: {e}")

        return predictions

    def run_single(
        self,
        params: np.ndarray,
        power: float,
    ) -> Tuple[float, float, float]:
        """
        运行单个功率点的仿真

        Parameters
        ----------
        params : np.ndarray
            参数（active 参数向量）
        power : float
            激光功率 (W)

        Returns
        -------
        tuple
            (width, depth, area) 单位: μm, μm, μm²
        """
        result = self.run(params, np.array([power]))
        return tuple(result[0])

    def evaluate_objective(
        self,
        params: np.ndarray,
        exp_data: np.ndarray,
    ) -> float:
        """
        评估目标函数（Bayes 统一 NRMSE%）

        Parameters
        ----------
        params : np.ndarray
            参数（active 参数向量）
        exp_data : np.ndarray, shape (n, 4)
            实验数据 [power, width, depth, area]

        Returns
        -------
        float
            NRMSE 百分比
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
        计算展平的残差向量（用于 scipy.optimize.least_squares）

        这里返回 Bayes 统一目标对应的归一化残差（支持 output_weights）。
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
        """重置评估计数器"""
        self.eval_count = 0
