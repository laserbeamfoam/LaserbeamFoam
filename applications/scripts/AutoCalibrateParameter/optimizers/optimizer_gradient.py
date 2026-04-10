"""
梯度优化器

使用 scipy.optimize.least_squares 进行梯度优化（点估计）
"""

from __future__ import annotations

import json
import threading
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, Optional

import numpy as np
from scipy.optimize import least_squares

from .base_optimizer import BaseOptimizer, OptimizationResult
from ..simulation.objective import compute_bayes_nrmse_percent, compute_bayes_normalized_residuals

if TYPE_CHECKING:
    from ..config.gradient_config import GradientConfig
    from ..simulation.simulation_runner import SimulationRunner


class GradientOptimizer(BaseOptimizer):
    """
    梯度优化器 (scipy.optimize.least_squares)

    特点：
    - 基于有限差分的梯度估计
    - 多起点并行优化
    - 返回点估计

    Attributes
    ----------
    config : GradientConfig
        梯度优化配置
    runner : SimulationRunner
        仿真运行器
    """

    def __init__(self, config: "GradientConfig", runner: "SimulationRunner"):
        """
        初始化梯度优化器

        Parameters
        ----------
        config : GradientConfig
            梯度优化配置
        runner : SimulationRunner
            仿真运行器
        """
        super().__init__(config)
        self.config = config
        self.runner = runner

    def get_name(self) -> str:
        return "Gradient"

    def optimize(self, exp_data: np.ndarray) -> OptimizationResult:
        """
        执行梯度优化

        Parameters
        ----------
        exp_data : np.ndarray, shape (n, 4)
            实验数据 [power, width, depth, area]

        Returns
        -------
        OptimizationResult
        """
        cfg = self.config

        print(f"\n{'='*60}")
        print(f"梯度优化 (scipy.optimize.least_squares)")
        print(f"{'='*60}")
        print(f"方向数: {cfg.n_directions}, 最大评估: {cfg.max_nfev}")
        print(f"初始化策略: {cfg.init_strategy}")
        print(f"并行 workers: {cfg.max_workers}")
        print()

        # 运行目录
        cfg.runs_root.mkdir(parents=True, exist_ok=True)

        # 存储所有方向的结果
        all_results: List[Dict] = []

        # 启动监控线程
        stop_event = threading.Event()
        monitor_thread = threading.Thread(
            target=self._monitor_loop,
            args=(stop_event,),
            daemon=True,
        )
        monitor_thread.start()

        try:
            # 依次运行每个方向（为简化，不使用进程池）
            for direction_id in range(1, cfg.n_directions + 1):
                print(f"\n方向 {direction_id}/{cfg.n_directions}")
                result = self._optimize_single_direction(direction_id, exp_data)
                all_results.append(result)
                print(f"  完成: NRMSE = {result['cost']:.4f}%")

        finally:
            stop_event.set()
            monitor_thread.join(timeout=2.0)

        # 选择最优结果
        if not all_results:
            return OptimizationResult(
                method="Gradient",
                best_params=np.zeros(len(self.runner.active_names)),
                best_cost=float("inf"),
                n_evaluations=0,
                history={},
                message="所有方向优化失败",
            )

        best = min(all_results, key=lambda r: r["cost"])
        total_evals = sum(r.get("nfev", 0) for r in all_results)

        print(f"\n最优方向: {best['direction_id']}, NRMSE = {best['cost']:.4f}%")

        return OptimizationResult(
            method="Gradient",
            best_params=np.array(best["params"]),
            best_cost=best["cost"],
            n_evaluations=total_evals,
            history={"all_results": all_results},
            message=f"梯度优化完成，最优方向 {best['direction_id']}",
            extra={"all_directions": all_results},
        )

    def _optimize_single_direction(
        self,
        direction_id: int,
        exp_data: np.ndarray,
    ) -> Dict:
        """
        单方向优化

        Parameters
        ----------
        direction_id : int
            方向 ID
        exp_data : np.ndarray
            实验数据

        Returns
        -------
        dict
            优化结果
        """
        cfg = self.config
        bounds = self.runner.get_active_bounds()

        # 创建工作目录
        work_dir = cfg.runs_root / f"direction_{direction_id:03d}"
        work_dir.mkdir(parents=True, exist_ok=True)

        # 初始点
        x0 = self._get_initial_point(direction_id, bounds)
        lb = np.array([b[0] for b in bounds])
        ub = np.array([b[1] for b in bounds])

        # 历史记录
        cost_history: List[float] = []
        params_history: List[List[float]] = []

        # 保存进度文件
        progress_file = work_dir / "progress.json"

        def residual_func(params: np.ndarray) -> np.ndarray:
            """残差函数（Bayes 统一归一化残差）"""
            power_points = exp_data[:, 0]
            observations = exp_data[:, 1:4]

            predictions = self.runner.run(params, power_points)

            weights = np.asarray(cfg.output_weights, dtype=float)
            residuals = compute_bayes_normalized_residuals(
                predictions,
                observations,
                weights=weights,
            ).flatten()

            cost, _ = compute_bayes_nrmse_percent(
                predictions,
                observations,
                weights=weights,
            )
            cost_history.append(cost)
            params_history.append(params.tolist())

            # 更新进度
            self._save_progress(
                progress_file,
                direction_id,
                len(cost_history),
                cost,
                min(cost_history),
            )

            if len(cost_history) % cfg.log_every_eval == 0:
                print(f"    [eval {len(cost_history)}] NRMSE = {cost:.4f}%")

            return residuals

        # 运行优化
        result = least_squares(
            residual_func,
            x0,
            bounds=(lb, ub),
            max_nfev=cfg.max_nfev,
            ftol=cfg.ftol,
            xtol=cfg.xtol,
            diff_step=cfg.diff_step,
        )

        final_predictions = self.runner.run(result.x, exp_data[:, 0])
        final_cost, _ = compute_bayes_nrmse_percent(
            final_predictions,
            exp_data[:, 1:4],
            weights=np.asarray(cfg.output_weights, dtype=float),
        )

        # 保存结果
        result_dict = {
            "direction_id": direction_id,
            "params": result.x.tolist(),
            "cost": final_cost,
            "nfev": result.nfev,
            "status": int(result.status),
            "message": result.message,
            "x0": x0.tolist(),
            "cost_history": cost_history,
        }

        result_file = work_dir / "result.json"
        with open(result_file, "w") as f:
            json.dump(result_dict, f, indent=2)

        return result_dict

    def _get_initial_point(
        self,
        direction_id: int,
        bounds: List,
    ) -> np.ndarray:
        """
        生成初始点

        Parameters
        ----------
        direction_id : int
            方向 ID
        bounds : list
            参数边界

        Returns
        -------
        np.ndarray
            初始参数
        """
        cfg = self.config
        n_params = len(bounds)

        if cfg.init_strategy == "legacy_single":
            # 所有方向使用中点
            return np.array([(b[0] + b[1]) / 2 for b in bounds])

        elif cfg.init_strategy == "legacy_multi":
            # 随机初始化
            rng = np.random.default_rng(seed=1000 + direction_id)
            return np.array([rng.uniform(b[0], b[1]) for b in bounds])

        elif cfg.init_strategy == "directional":
            # 从中点出发，沿随机方向偏移
            center = np.array([(b[0] + b[1]) / 2 for b in bounds])
            scales = np.array([(b[1] - b[0]) / 2 for b in bounds])

            rng = np.random.default_rng(seed=1000 + direction_id)
            direction = rng.normal(0, 1, size=n_params)
            direction = direction / (np.linalg.norm(direction) + 1e-10)

            x0 = center + cfg.direction_scale * scales * direction

            # 裁剪到边界内
            for i in range(n_params):
                x0[i] = np.clip(x0[i], bounds[i][0], bounds[i][1])

            return x0

        else:
            raise ValueError(f"未知初始化策略: {cfg.init_strategy}")

    def _save_progress(
        self,
        filepath: Path,
        direction_id: int,
        n_eval: int,
        current_cost: float,
        best_cost: float,
    ) -> None:
        """保存进度信息"""
        data = {
            "direction_id": direction_id,
            "n_eval": n_eval,
            "current_cost": current_cost,
            "best_cost": best_cost,
            "ts": time.time(),
        }
        with open(filepath, "w") as f:
            json.dump(data, f)

    def _monitor_loop(self, stop_event: threading.Event) -> None:
        """监控线程"""
        cfg = self.config
        poll_interval = max(0.1, cfg.monitor_interval_sec)

        while not stop_event.is_set():
            # 收集所有方向的进度
            rows = []
            for i in range(1, cfg.n_directions + 1):
                progress_file = cfg.runs_root / f"direction_{i:03d}" / "progress.json"
                if progress_file.exists():
                    try:
                        with open(progress_file) as f:
                            d = json.load(f)
                        rows.append((
                            i,
                            d.get("n_eval", 0),
                            d.get("current_cost"),
                            d.get("best_cost"),
                        ))
                    except Exception:
                        pass

            # 这里可以打印汇总信息（目前省略以减少输出）

            time.sleep(poll_interval)
