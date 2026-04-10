"""
贝叶斯优化器

使用 scikit-optimize 进行贝叶斯优化（点估计）
"""

from __future__ import annotations

import json
import os
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Tuple

import numpy as np

try:
    from skopt import Optimizer
    from skopt.space import Real
except ImportError:
    Optimizer = None
    Real = None

try:
    from pyDOE import lhs
except ImportError:
    lhs = None

import pandas as pd

from .base_optimizer import BaseOptimizer, OptimizationResult
from ..simulation.objective import compute_bayes_nrmse_percent
from ..utils.param_registry import (
    get_param_csv_candidates,
    get_param_csv_column,
    get_param_csv_column_map,
)

if TYPE_CHECKING:
    from ..config.bayes_config import BayesConfig
    from ..simulation.simulation_runner import SimulationRunner


def _evaluate_single_point(args: Tuple) -> Tuple[int, List[float], float]:
    """
    单点评估函数（用于并行）

    Parameters
    ----------
    args : tuple
        (job_id, params, config_dict, exp_data, output_weights)

    Returns
    -------
    tuple
        (job_id, params, cost)
    """
    job_id, params, config_dict, exp_data, output_weights = args

    # 重建配置和运行器
    from ..config.bayes_config import BayesConfig
    from ..simulation.simulation_runner import SimulationRunner
    from ..simulation.objective import compute_bayes_nrmse_percent

    config = BayesConfig(**config_dict)
    runner = SimulationRunner(config)

    # 评估
    power_points = exp_data[:, 0]
    observations = exp_data[:, 1:4]

    predictions = runner.run(np.array(params), power_points, job_id=job_id)

    cost, _ = compute_bayes_nrmse_percent(
        predictions,
        observations,
        weights=np.asarray(output_weights, dtype=float),
    )

    return job_id, params, cost


class BayesianOptimizer(BaseOptimizer):
    """
    贝叶斯优化器 (scikit-optimize)

    特点：
    - 高斯过程代理 + 采集函数
    - 返回点估计（最优解）
    - 支持批量并行评估

    Attributes
    ----------
    config : BayesConfig
        贝叶斯优化配置
    runner : SimulationRunner
        仿真运行器
    """

    def __init__(self, config: "BayesConfig", runner: "SimulationRunner"):
        """
        初始化贝叶斯优化器

        Parameters
        ----------
        config : BayesConfig
            贝叶斯优化配置
        runner : SimulationRunner
            仿真运行器
        """
        if Optimizer is None:
            raise ImportError("需要安装 scikit-optimize: pip install scikit-optimize")

        super().__init__(config)
        self.config = config
        self.runner = runner

        # 输出权重
        self.output_weights = np.array(config.output_weights)

    def get_name(self) -> str:
        return "Bayesian"

    def _get_history_csv_path(self) -> Path:
        """获取唯一的 history CSV 路径（统一放在 runs_root 下）。"""
        history_name = "bayes_history.csv"
        if getattr(self.config, "warmstart_csv", None):
            configured_raw = str(self.config.warmstart_csv)
            configured_name = Path(configured_raw).name
            if configured_raw != configured_name:
                print(f"警告: warmstart_csv='{configured_raw}' 包含路径，已规范为文件名 '{configured_name}' 并存储到 runs_root")
            if configured_name:
                history_name = configured_name
        return self.config.runs_root / history_name

    def _get_progress_file_path(self, job_id: int) -> Path:
        """获取指定 job 的进度文件路径（统一放在 runs_root 下）。"""
        return self.config.runs_root / f"bayes_progress_job{job_id}.csv"

    def _migrate_legacy_files(self) -> None:
        """将历史遗留的 case 根目录文件迁移到 runs_root，保证只保留 runs 一份。"""
        self.config.runs_root.mkdir(parents=True, exist_ok=True)

        # 迁移 bayes_history.csv（若 runs 中不存在）
        history_in_runs = self._get_history_csv_path()
        legacy_history = self.config.case_dir / history_in_runs.name
        if legacy_history.exists() and legacy_history != history_in_runs:
            if history_in_runs.exists():
                try:
                    legacy_history.unlink()
                    print(f"清理旧路径 history 文件: {legacy_history}")
                except Exception as e:
                    print(f"警告: 检测到旧路径 history {legacy_history}，无法删除({e})，将继续使用 {history_in_runs}")
            else:
                shutil.move(str(legacy_history), str(history_in_runs))
                print(f"迁移历史文件到 runs: {history_in_runs}")

        # 迁移 bayes_progress_job*.csv
        for legacy_progress in sorted(self.config.case_dir.glob("bayes_progress_job*.csv")):
            target_progress = self.config.runs_root / legacy_progress.name
            if legacy_progress == target_progress:
                continue
            if target_progress.exists():
                try:
                    legacy_progress.unlink()
                    print(f"清理旧进度文件: {legacy_progress}")
                except Exception as e:
                    print(f"警告: 无法清理旧进度文件 {legacy_progress}: {e}")
            else:
                shutil.move(str(legacy_progress), str(target_progress))
                print(f"迁移进度文件到 runs: {target_progress}")

    def optimize(self, exp_data: np.ndarray) -> OptimizationResult:
        """
        执行贝叶斯优化

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
        print(f"贝叶斯优化 (scikit-optimize)")
        print(f"{'='*60}")
        print(f"采集函数: {cfg.acq_func}")
        print(f"初始点: {cfg.n_initial_points}, 批次: {cfg.n_batches}, 批大小: {cfg.batch_size}")
        print()

        # 定义搜索空间（使用runner的active参数边界）
        bounds = self.runner.get_active_bounds()
        param_names = self.runner.active_names
        dimensions = [
            Real(low, high, name=name)
            for (low, high), name in zip(bounds, param_names)
        ]

        # 初始化 skopt 优化器
        optimizer = Optimizer(
            dimensions=dimensions,
            base_estimator="GP",
            acq_func=cfg.acq_func,
            acq_func_kwargs={"xi": cfg.xi, "kappa": cfg.kappa},
            n_initial_points=0,  # 我们自己做初始采样
            random_state=42,
        )

        all_x: List[List[float]] = []
        all_y: List[float] = []
        all_predictions: List[Optional[np.ndarray]] = []  # 保存所有预测结果
        best_cost = float("inf")

        # 获取功率点用于保存 CSV
        power_points = exp_data[:, 0]

        # 统一迁移历史遗留文件到 runs_root
        self._migrate_legacy_files()

        # 检查是否从 runs 下的 CSV 热启动
        history_csv_path = self._get_history_csv_path()
        if history_csv_path.exists():
            try:
                # Load history and predictions
                warmstart_x, warmstart_y, warmstart_preds = self._load_history_csv(str(history_csv_path))
                if warmstart_x and warmstart_y:
                    all_x.extend(warmstart_x)
                    all_y.extend(warmstart_y)
                    # Use loaded predictions instead of None
                    all_predictions.extend(warmstart_preds)
                    best_cost = min(all_y)
                    print(f"从 CSV 热启动: {len(warmstart_y)} 个样本, 最优 RMSE = {best_cost:.2f}")
            except Exception as e:
                print(f"加载 CSV 热启动数据失败: {e}")

        # 检查是否从 history.json 恢复 (Deprecated or secondary)
        # 检查是否从 history.json 恢复 (已弃用，只使用 CSV)
        # history_file = cfg.runs_root / "history.json"
        # if cfg.resume and history_file.exists():
        #     try:
        #         with open(history_file, "r") as f:
        #             history_data = json.load(f)
        #         resume_x = history_data.get("x_iters", [])
        #         resume_y = history_data.get("y_iters", [])
        #         if resume_x and resume_y:
        #             # 合并数据（避免重复）
        #             existing_set = set(tuple(x) for x in all_x)
        #             for x, y in zip(resume_x, resume_y):
        #                 if tuple(x) not in existing_set:
        #                     all_x.append(x)
        #                     all_y.append(y)
        #                     all_predictions.append(None)  # 恢复的数据没有预测结果
        #                     existing_set.add(tuple(x))
        #             if all_y:
        #                 best_cost = min(all_y)
        #             print(f"从 history.json 恢复: 当前共 {len(all_y)} 个样本, 最优 RMSE = {best_cost:.2f}")
        #     except Exception as e:
        #         print(f"恢复历史失败: {e}")

        # 向优化器提供历史数据
        if all_x and all_y:
            # Important: Ensure NO NaNs in all_y before telling optimizer
            valid_indices = [i for i, y in enumerate(all_y) if not np.isnan(y)]
            if len(valid_indices) < len(all_y):
                 print(f"警告: 忽略 {len(all_y) - len(valid_indices)} 个 NaN 目标值记录")
            
            clean_x = [all_x[i] for i in valid_indices]
            clean_y = [all_y[i] for i in valid_indices]

            # Clip historical points to current bounds (handles cases where bounds were widened/narrowed)
            bounds_list = [space_dim.bounds for space_dim in optimizer.space.dimensions]
            clipped_x = []
            n_clipped = 0
            for pt in clean_x:
                clipped_pt = [float(np.clip(v, lo, hi)) for v, (lo, hi) in zip(pt, bounds_list)]
                if clipped_pt != list(pt):
                    n_clipped += 1
                clipped_x.append(clipped_pt)
            if n_clipped > 0:
                print(f"警告: {n_clipped} 个历史点超出新边界，已裁剪至边界")

            if clipped_x and clean_y:
                optimizer.tell(clipped_x, clean_y)
                print(f"总历史数据: {len(clean_y)} 个有效样本")

        # 阶段 1: 初始采样
        n_initial_needed = max(0, cfg.n_initial_points - len(all_y))

        # 检测未完成的进度文件，优先恢复
        incomplete_jobs = self._detect_incomplete_jobs()
        if incomplete_jobs:
            print(f"\n检测到 {len(incomplete_jobs)} 个未完成的 jobs，优先继续评估")
            for job_info in incomplete_jobs:
                job_id = job_info['job_id']
                params = job_info['params']
                print(f"  [Job {job_id}] 继续评估: {params}")

                # 直接评估这个未完成的job
                n_before = len(all_y)
                batch_results = self._evaluate_batch(
                    [params], exp_data, n_before,
                    all_x=all_x, all_y=all_y, all_predictions=all_predictions,
                    power_points=exp_data[:, 0],
                    explicit_job_ids=[job_id]
                )

                # 更新历史
                for param_vals, cost, pred in batch_results:
                    all_x.append(param_vals)
                    all_y.append(cost)
                    all_predictions.append(pred)
                    if not np.isnan(cost) and (best_cost is None or cost < best_cost):
                        best_cost = cost

                # 减少需要采样的数量
                n_initial_needed = max(0, n_initial_needed - 1)

        if n_initial_needed > 0:
            print(f"\n阶段 1: 初始采样 ({n_initial_needed} 点)")
            initial_samples = self._initial_sampling(n_initial_needed, bounds)

            for batch_start in range(0, n_initial_needed, cfg.batch_size):
                batch_end = min(batch_start + cfg.batch_size, n_initial_needed)
                batch_samples = initial_samples[batch_start:batch_end]

                # 记录批次开始前的数量，用于获取新增结果
                n_before = len(all_y)

                # 调用 _evaluate_batch，它会实时保存每个模拟结果
                batch_results = self._evaluate_batch(
                    batch_samples, exp_data, n_before,
                    all_x=all_x, all_y=all_y, all_predictions=all_predictions,
                    power_points=power_points
                )

                # 更新最优值并记录日志
                for i, (params, cost, preds) in enumerate(batch_results):
                    if cost < best_cost:
                        best_cost = cost
                    self._log_iteration(np.array(params), cost, n_before + i + 1)

                print(f"  批次完成, 当前最优 RMSE = {best_cost:.2f}")
                self._generate_interim_plots(
                    all_x, all_y, all_predictions, exp_data, power_points,
                    n_initial=cfg.n_initial_points,
                )

        # 阶段 2: 贝叶斯优化
        print(f"\n阶段 2: 贝叶斯优化 ({cfg.n_batches} 批次)")

        # Re-tell only valid data just in case
        # (Though we did it before, newly added initial samples are automatically handled by loop or manual tell?)
        # skopt optimizer handles tell() incrementally. We don't need to re-tell everything.
        # But we DO need to tell the INITIAL SAMPLES we just ran if they weren't told.
        # However, we only told the `warmstart` data. The new `initial_samples` need to be told.
        # Since we use `ask()` in loop below, we must be careful.
        # Actually, standard usage is: tell initial points, then ask-tell loop.
        
        # Tell the initial samples we just ran
        if n_initial_needed > 0:
             # Identify the indices of the new samples in all_x/all_y
             # They are at the end: from len(all_y)-n_initial_needed to end?
             # No, we might have had failures.
             # Simply: we haven't told the optimizer about the new samples yet.
             # In skopt, `tell` returns a result object, but we manage state ourselves too.
             
             # We need to tell the optimizer about the samples obtained in Phase 1
             # We can just tell everything again? No, that duplicates.
             # We should tell the NEW samples.
             
             # The warmstart data was told. Now we tell the Phase 1 data.
             # Phase 1 added samples to all_x, all_y.
             # So we tell the slice that wasn't told.
             start_idx = len(all_x) - n_initial_needed # Assuming all succeeded?
             # That assumption is risky. Let's rely on checking what's new.
             # Actually, simpler: construct a clean list of EVERYTHING and re-initialize optimizer? 
             # No, optimizer has internal state (GP fit).
             
             # Better: just tell the new points from Phase 1.
             # Since we know exactly which ones they are (from batch_results).
             # We didn't accumulate batch_results from the loop easily.
             
             # Retrospective fix: We iterate batches in Phase 1. 
             # We should collect them and tell.
             pass 
             
             # Actually, let's just make sure we tell ALL valid points that haven't been told.
             # But optimizer.tell() is stateful.
             # Let's clean this up:
             # 1. New Optimizer is created.
             # 2. We tell warmstart data.
             # 3. We run Phase 1 (Initial Sampling).
             # 4. We MUST tell Phase 1 data to optimizer before Phase 2.
        
        # Let's grab the data generated in Phase 1
        # It's all_x/all_y minus the warmstart count.
        # But wait, we don't know the exact count easily unless we tracked it.
        # Let's count how many we told initially.
        n_told = 0
        if all_x and all_y:
             # We told everything available at start
             # But we filtered NaNs.
             pass 
        
        # To be safe and simple: 
        # For the Phase 1 points, we should iterate and tell them.
        # But doing it inside the loop is messy.
        
        # Correct approach:
        # We collected all_x/all_y. The optimizer only knows about the warmstart part.
        # We should tell the rest.
        
        # Let's re-fit the optimizer with ALL current data to be sure.
        # But skopt doesn't have a 'reset' easily.
        # Actually, we can just tell the new points.
        
        # IMPORTANT: The original code didn't tell Phase 1 points to optimizer explicitly!
        # It relied on `optimizer` being used in Phase 2? No.
        # If we don't tell Phase 1 points, Phase 2 starts with only warmstart info.
        # That's a bug in original code too, potentially.
        # Wait, usually `n_initial_points` in `Optimizer` handles this?
        # We set `n_initial_points=0` and do it manually.
        
        # So we MUST tell the points from Phase 1.
        # Let's identify them.
        # We can just check `optimizer.Xi` to see what it has (if accessible) or just track manually.
        
        # Simplest hack: We know `all_x` contains everything.
        # We can create a NEW optimizer for Phase 2, fed with ALL history.
        optimizer = Optimizer(
            dimensions=dimensions,
            base_estimator="GP",
            acq_func=cfg.acq_func,
            acq_func_kwargs={"xi": cfg.xi, "kappa": cfg.kappa},
            n_initial_points=0,
            random_state=42,
        )
        
        # Tell ALL valid history
        valid_indices = [i for i, y in enumerate(all_y) if not np.isnan(y)]
        clean_x = [all_x[i] for i in valid_indices]
        clean_y = [all_y[i] for i in valid_indices]

        # Clip to current bounds (handles cases where bounds were widened/narrowed)
        bounds_list = [d.bounds for d in optimizer.space.dimensions]
        clean_x = [[float(np.clip(v, lo, hi)) for v, (lo, hi) in zip(pt, bounds_list)] for pt in clean_x]

        if clean_x and clean_y:
            optimizer.tell(clean_x, clean_y)


        for batch_idx in range(cfg.n_batches):
            # 获取候选点
            candidates = optimizer.ask(n_points=cfg.batch_size)
            candidates = [list(c) for c in candidates]

            print(f"\n批次 {batch_idx + 1}/{cfg.n_batches}")

            # 记录批次开始前的数量
            n_before = len(all_y)

            # 调用 _evaluate_batch，它会实时保存每个模拟结果
            batch_results = self._evaluate_batch(
                candidates, exp_data, n_before,
                all_x=all_x, all_y=all_y, all_predictions=all_predictions,
                power_points=power_points
            )

            # 收集本批次的结果用于告知优化器
            batch_x = [r[0] for r in batch_results]
            batch_y = [r[1] for r in batch_results]

            # 更新最优值并记录日志
            for i, (params, cost, preds) in enumerate(batch_results):
                if cost < best_cost:
                    best_cost = cost
                self._log_iteration(np.array(params), cost, n_before + i + 1)

            # Only tell valid results
            valid_batch_indices = [i for i, y in enumerate(batch_y) if not np.isnan(y)]
            if valid_batch_indices:
                clean_batch_x = [batch_x[i] for i in valid_batch_indices]
                clean_batch_y = [batch_y[i] for i in valid_batch_indices]
                optimizer.tell(clean_batch_x, clean_batch_y)
            
            if len(valid_batch_indices) < len(batch_y):
                 print(f"  警告: 本批次忽略 {len(batch_y) - len(valid_batch_indices)} 个失败/NaN 结果")

            print(f"  批次 {batch_idx + 1} 完成, 最优 RMSE = {best_cost:.2f}")
            self._generate_interim_plots(
                all_x, all_y, all_predictions, exp_data, power_points,
                n_initial=cfg.n_initial_points,
            )

        # 提取最优解
        best_idx = int(np.argmin(all_y))
        best_params = np.array(all_x[best_idx])

        return OptimizationResult(
            method="Bayesian",
            best_params=best_params,
            best_cost=all_y[best_idx],
            n_evaluations=len(all_y),
            history={"params": all_x, "costs": all_y},
            message=f"贝叶斯优化完成，共 {len(all_y)} 次评估",
        )

    def _generate_interim_plots(
        self,
        all_x: List[List[float]],
        all_y: List[float],
        all_predictions: List,
        exp_data: np.ndarray,
        power_points: np.ndarray,
        n_initial: int = 0,
    ) -> None:
        """每批次完成后生成中间可视化图（覆盖更新，不累积文件）"""
        from ..postprocess.comparison import save_comparison_plot, plot_predictions_vs_experiments
        from ..postprocess.convergence import save_convergence_plot

        valid_indices = [i for i, y in enumerate(all_y) if not np.isnan(y)]
        if not valid_indices:
            return

        best_idx = min(valid_indices, key=lambda i: all_y[i])
        best_cost = all_y[best_idx]
        best_preds = all_predictions[best_idx]

        if best_preds is None:
            return

        # 只取 width/depth 两列，忽略 area
        observations = exp_data[:, 1:3]
        best_preds_wd = best_preds[:, :2]

        plots_dir = self.config.runs_root / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)

        n_evals = len(valid_indices)

        try:
            save_comparison_plot(
                plots_dir,
                power_points,
                observations,
                best_preds_wd,
                output_names=["Width", "Depth"],
                output_units=["um", "um"],
                nrmse=best_cost,
                title=f"Best result (n={n_evals}, NRMSE={best_cost:.2f}%)",
            )
            plot_predictions_vs_experiments(
                plots_dir,
                observations,
                best_preds_wd,
                output_names=["Width", "Depth"],
                output_units=["um", "um"],
                title=f"Predictions vs Experiments (n={n_evals})",
            )
        except Exception as e:
            print(f"  [警告] 生成对比图失败: {e}")

        try:
            valid_costs = [all_y[i] for i in valid_indices]
            save_convergence_plot(
                plots_dir, valid_costs, n_initial=n_initial,
                title="Optimization Convergence",
            )
        except Exception as e:
            print(f"  [警告] 生成收敛图失败: {e}")

        try:
            from ..postprocess.convergence import plot_parameter_evolution
            valid_params = [all_x[i] for i in valid_indices]
            plot_parameter_evolution(
                plots_dir,
                [np.array(p) for p in valid_params],
                self.runner.active_names,
                title="Parameter Evolution",
            )
        except Exception as e:
            print(f"  [警告] 生成参数演化图失败: {e}")

        print(f"  [图表] 已更新: {plots_dir}")

    def _initial_sampling(
        self,
        n_points: int,
        bounds: List[Tuple[float, float]],
    ) -> List[List[float]]:
        """初始 LHS 采样"""
        rng_seed = 42 if self.config.seed is None else int(self.config.seed)

        if lhs is None:
            # 后备：随机采样
            np.random.seed(rng_seed)
            samples = []
            for _ in range(n_points):
                sample = [
                    np.random.uniform(low, high)
                    for low, high in bounds
                ]
                samples.append(sample)
            return samples

        np.random.seed(rng_seed)
        lhs_samples = lhs(len(bounds), n_points)
        samples = []
        for row in lhs_samples:
            sample = [
                bounds[i][0] + row[i] * (bounds[i][1] - bounds[i][0])
                for i in range(len(bounds))
            ]
            samples.append(sample)
        return samples

    def _evaluate_batch(
        self,
        candidates: List[List[float]],
        exp_data: np.ndarray,
        offset: int,
        all_x: List[List[float]] = None,
        all_y: List[float] = None,
        all_predictions: List = None,
        power_points: np.ndarray = None,
        explicit_job_ids: List[int] = None,
    ) -> List[Tuple[List[float], float, Optional[np.ndarray]]]:
        """
        评估一批候选点
        """
        results = []

        for idx, params in enumerate(candidates):
            if explicit_job_ids:
                job_id = explicit_job_ids[idx]
            else:
                job_id = offset + idx + 1

            # 检查是否存在部分完成的进度文件
            progress_file = self._get_progress_file_path(job_id)
            partial_preds = None
            start_power_idx = 0

            if progress_file.exists():
                try:
                    partial_preds, start_power_idx = self._load_partial_progress(
                        job_id, params, power_points
                    )
                    if partial_preds is not None and start_power_idx > 0:
                        print(f"  [Job {job_id}] 检测到部分进度，从功率点 {start_power_idx+1}/{len(power_points)} 继续")
                except Exception as e:
                    print(f"  [Job {job_id}] 无法加载部分进度: {e}，从头开始")
                    partial_preds = None
                    start_power_idx = 0
            # 动态构建参数日志字符串
            param_strs = []
            for i, name in enumerate(self.runner.active_names):
                if name == "sigma":
                    param_strs.append(f"sigma={params[i]:.4f}")
                elif name == "marangoni":
                    param_strs.append(f"marangoni={params[i]:.2e}")
                elif name == "substrate_temp":
                    param_strs.append(f"T_s={params[i]:.1f}K")
                elif name == "absorptivity":
                    param_strs.append(f"absorptivity={params[i]:.3f}")
                elif name == "recoilCoeff":
                    param_strs.append(f"recoilCoeff={params[i]:.3f}")
                elif name == "radius_flavour":
                    param_strs.append(f"Radius_Flavour={params[i]:.3f}")
                else:
                    param_strs.append(f"{name}={params[i]:.4g}")
            print(f"  [Job {job_id}] 参数: {', '.join(param_strs)}")

            try:
                exp_power_points = exp_data[:, 0]
                observations = exp_data[:, 1:4]

                # 创建回调函数，在每个功率点完成后保存中间结果
                partial_predictions_storage = {'data': None}

                def on_power_complete(power_idx: int, power_val: float, preds_so_far: np.ndarray):
                    """每个功率点完成后的回调"""
                    partial_predictions_storage['data'] = preds_so_far.copy()

                    # 计算当前已完成功率点的部分 NRMSE（与 Bayes 目标保持一致）
                    n_completed = power_idx + 1
                    partial_obs = observations[:n_completed, :]
                    partial_preds = preds_so_far

                    partial_cost, _ = compute_bayes_nrmse_percent(
                        partial_preds,
                        partial_obs,
                        weights=np.asarray(self.output_weights, dtype=float),
                    )

                    print(f"  [Job {job_id}] 功率点 {power_idx+1}/{len(exp_power_points)} 完成 "
                          f"(P={power_val:.0f}W), 部分NRMSE={partial_cost:.2f}%")

                    # 实时保存中间结果到临时CSV
                    if all_x is not None and power_points is not None:
                        self._save_partial_progress(
                            job_id, params, partial_preds, power_points[:n_completed],
                            partial_cost, n_completed, len(exp_power_points)
                        )

                # 如果有部分完成的数据，从断点继续
                if partial_preds is not None and start_power_idx > 0:
                    # 初始化预测数组，填入已完成的数据
                    predictions = partial_preds.copy()

                    # 只运行剩余的功率点
                    remaining_powers = exp_power_points[start_power_idx:]
                    if len(remaining_powers) > 0:
                        # 创建增强的回调函数，传递完整的预测数组
                        def enhanced_callback(power_idx: int, power_val: float, partial_remaining: np.ndarray):
                            """包装回调函数，传递完整数据"""
                            # 更新完整预测数组
                            predictions[power_idx, :] = partial_remaining[-1, :]

                            # 调用原始回调，传递完整数据
                            if on_power_complete is not None:
                                on_power_complete(power_idx, power_val, predictions[:power_idx+1].copy())

                        remaining_preds = self.runner.run(
                            np.array(params), remaining_powers, job_id=job_id,
                            on_power_complete=enhanced_callback if on_power_complete else None,
                            start_index=start_power_idx  # 传递起始索引
                        )
                        # 填充剩余的预测结果
                        predictions[start_power_idx:] = remaining_preds
                else:
                    # 从头开始运行所有功率点
                    predictions = self.runner.run(
                        np.array(params), exp_power_points, job_id=job_id,
                        on_power_complete=on_power_complete
                    )

                cost, nrmse_per_output = compute_bayes_nrmse_percent(
                    predictions,
                    observations,
                    weights=np.asarray(self.output_weights, dtype=float),
                )
                output_names = ['width', 'depth', 'area']
                weights = np.asarray(self.output_weights, dtype=float)
                nrmse_str = ', '.join([
                    f"{name}={nrmse_per_output[i]:.2f}%{'*' if weights[i] > 0 else ''}"
                    for i, name in enumerate(output_names)
                ])
                print(f"  [Job {job_id}] NRMSE = {cost:.2f}% ({nrmse_str})")
                results.append((params, cost, predictions))

            except Exception as e:
                print(f"  [Job {job_id}] 失败: {e}")
                results.append((params, np.nan, None)) # Use np.nan for failure instead of 1e16 to handle properly

            # 实时保存
            if all_x is not None and all_y is not None and all_predictions is not None:
                all_x.append(results[-1][0])
                all_y.append(results[-1][1])
                all_predictions.append(results[-1][2])
                self._save_history(all_x, all_y)
                if power_points is not None:
                    self._save_detailed_csv(all_x, all_y, all_predictions, power_points)
                print(f"  [Job {job_id}] 数据已保存 (共 {len(all_y)} 条记录)")

                # 删除临时进度文件（如果存在）
                progress_file = self._get_progress_file_path(job_id)
                if progress_file.exists():
                    try:
                        progress_file.unlink()
                        print(f"  [Job {job_id}] 临时进度文件已删除")
                    except Exception as e:
                        print(f"  [警告] 无法删除临时进度文件: {e}")

        return results

    def _save_history(
        self,
        all_x: List[List[float]],
        all_y: List[float],
    ) -> None:
        """保存优化历史 (已弃用，不再保存 json)"""
        pass
        # history_file = self.config.runs_root / "history.json"
        # ... (removed)

    def _save_detailed_csv(
        self,
        all_x: List[List[float]],
        all_y: List[float],
        all_predictions: List[Optional[np.ndarray]],
        power_points: np.ndarray,
    ) -> None:
        """
        保存详细的优化历史到 CSV 文件（仅 runs_root 下唯一文件）
        """
        rows = []
        for i, (params, cost, preds) in enumerate(zip(all_x, all_y, all_predictions)):
            # 动态构建参数列（根据active_names）
            row = {}
            for j, name in enumerate(self.runner.active_names):
                row[get_param_csv_column(name)] = params[j]

            # 添加固定参数（用于记录完整配置）
            for name, value in self.runner.fixed_dict.items():
                row[get_param_csv_column(name)] = value

            if preds is not None:
                for j, power in enumerate(power_points):
                    row[f"width_um_{int(power)}W"] = preds[j, 0]
                    row[f"depth_um_{int(power)}W"] = preds[j, 1]
                    row[f"area_um2_{int(power)}W"] = preds[j, 2]
            else:
                for power in power_points:
                    row[f"width_um_{int(power)}W"] = np.nan
                    row[f"depth_um_{int(power)}W"] = np.nan
                    row[f"area_um2_{int(power)}W"] = np.nan

            row["objective"] = cost 
            row["status"] = "ok" if (cost is not None and not np.isnan(cost) and cost < 1e6) else "failed"

            rows.append(row)

        df = pd.DataFrame(rows)

        csv_file = self._get_history_csv_path()
        csv_file.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(csv_file, index=False)

    @staticmethod
    def _resolve_param_csv_column(df: pd.DataFrame, param_name: str) -> str:
        """解析参数列名（优先新列名，兼容旧列名）。"""
        for candidate in get_param_csv_candidates(param_name):
            if candidate in df.columns:
                return candidate
        raise ValueError(f"CSV 缺少参数列: {param_name}, 候选列={get_param_csv_candidates(param_name)}")

    def _save_partial_progress(
        self,
        job_id: int,
        params: List[float],
        predictions: np.ndarray,
        power_points: np.ndarray,
        partial_cost: float,
        n_completed: int,
        n_total: int,
    ) -> None:
        """
        保存部分完成的进度到临时文件

        Parameters
        ----------
        job_id : int
            任务ID
        params : List[float]
            参数列表
        predictions : np.ndarray
            已完成的预测结果
        power_points : np.ndarray
            已完成的功率点
        partial_cost : float
            部分NRMSE
        n_completed : int
            已完成的功率点数量
        n_total : int
            总功率点数量
        """
        # 保存到临时进度文件
        progress_file = self._get_progress_file_path(job_id)

        # 动态构建参数行（与_save_detailed_csv保持一致）
        row = {"job_id": job_id}
        for j, name in enumerate(self.runner.active_names):
            row[get_param_csv_column(name)] = params[j]

        # 添加固定参数
        for name, value in self.runner.fixed_dict.items():
            row[get_param_csv_column(name)] = value

        row["progress"] = f"{n_completed}/{n_total}"
        row["partial_nrmse"] = partial_cost
        row["timestamp"] = pd.Timestamp.now().strftime("%Y/%m/%d %H:%M")

        # 添加已完成的功率点数据
        for i, power in enumerate(power_points):
            row[f"width_um_{int(power)}W"] = predictions[i, 0]
            row[f"depth_um_{int(power)}W"] = predictions[i, 1]
            row[f"area_um2_{int(power)}W"] = predictions[i, 2]

        df = pd.DataFrame([row])
        df.to_csv(progress_file, index=False)

        print(f"  [Job {job_id}] 进度已保存到 {progress_file.name}")

    def _detect_incomplete_jobs(self) -> List[Dict]:
        """
        检测未完成的 jobs（有 bayes_progress_job*.csv 但不在 bayes_history.csv 中）

        Returns
        -------
        list of dict
            每个字典包含 {'job_id': int, 'params': list, 'progress': str}
        """
        incomplete_jobs = []

        # 加载实验数据，获取总功率点数
        from ..data import load_experiment_data
        try:
            exp_data = load_experiment_data(self.config.exp_csv)
            total_power_points = len(exp_data)
            power_points = exp_data[:, 0]
        except Exception as e:
            print(f"  [警告] 无法加载实验数据: {e}")
            return incomplete_jobs

        # 扫描 runs_root 下的进度文件
        self.config.runs_root.mkdir(parents=True, exist_ok=True)
        for progress_file in sorted(self.config.runs_root.glob("bayes_progress_job*.csv")):
            try:
                df = pd.read_csv(progress_file)
                if len(df) == 0:
                    continue

                row = df.iloc[0]
                job_id = int(row['job_id'])

                # 提取参数（按 active_names 顺序）
                params = []
                for param_name in self.runner.active_names:
                    csv_col = self._resolve_param_csv_column(df, param_name)
                    if csv_col in row:
                        params.append(float(row[csv_col]))
                    else:
                        raise ValueError(f"进度文件缺少参数: {param_name}")

                # 通过检查实际的功率点数据来判断进度（不依赖progress字段）
                completed = 0
                for i, power in enumerate(power_points):
                    w_col = f"width_um_{int(power)}W"
                    d_col = f"depth_um_{int(power)}W"
                    a_col = f"area_um2_{int(power)}W"

                    # 检查这个功率点的数据是否存在且有效
                    if all(col in row for col in [w_col, d_col, a_col]):
                        if not any(pd.isna(row[col]) for col in [w_col, d_col, a_col]):
                            completed = i + 1
                        else:
                            break  # 遇到第一个无效数据就停止
                    else:
                        break  # 遇到第一个缺失列就停止

                # 只恢复部分完成的 jobs（0 < completed < total_power_points）
                if 0 < completed < total_power_points:
                    print(f"  [发现未完成Job] Job {job_id}: 已完成 {completed}/{total_power_points} 个功率点")
                    incomplete_jobs.append({
                        'job_id': job_id,
                        'params': params,
                        'progress': f"{completed}/{total_power_points}",
                        'progress_file': progress_file
                    })

            except Exception as e:
                print(f"  [警告] 读取进度文件 {progress_file.name} 失败: {e}")
                continue

        return incomplete_jobs

    def _load_partial_progress(
        self,
        job_id: int,
        params: List[float],
        power_points: np.ndarray,
    ) -> Tuple[Optional[np.ndarray], int]:
        """
        加载部分完成的进度

        Parameters
        ----------
        job_id : int
            任务ID
        params : List[float]
            当前参数（用于验证）
        power_points : np.ndarray
            所有功率点

        Returns
        -------
        tuple
            (partial_predictions, start_index)
            - partial_predictions: 已完成的预测数据，形状 (n_powers, 3)
            - start_index: 下一个需要运行的功率点索引
        """
        import pandas as pd

        progress_file = self._get_progress_file_path(job_id)
        if not progress_file.exists():
            return None, 0

        df = pd.read_csv(progress_file)
        if len(df) == 0:
            return None, 0

        row = df.iloc[0]

        # 验证参数是否匹配（可选，防止恢复错误的job）
        param_match = True
        for i, name in enumerate(self.runner.active_names):
            csv_name = self._resolve_param_csv_column(df, name)
            if csv_name in row and abs(row[csv_name] - params[i]) > 1e-6:
                param_match = False
                break

        if not param_match:
            print(f"  [警告] Job {job_id} 参数不匹配，忽略进度文件")
            return None, 0

        # 解析 progress 字段 (格式: "1/3" 表示已完成1个，共3个)
        progress_str = str(row.get("progress", "0/0"))
        try:
            # 尝试标准格式 "2/3"
            completed, total = map(int, progress_str.split("/"))
        except:
            # 尝试处理 Excel 日期格式 "3-Jan" -> 1/3 (1月3日 -> 月份=1, 日期=3)
            try:
                import re
                # 匹配 "数字-月份" 格式
                match = re.match(r'(\d+)-([A-Za-z]+)', progress_str)
                if match:
                    day = int(match.group(1))  # 日期（总数）
                    month_str = match.group(2).capitalize()  # 月份名称

                    # 月份映射：Excel 将 1/3 转为 1月3日 (3-Jan)
                    month_map = {
                        "Jan": 1, "Feb": 2, "Mar": 3, "Apr": 4, "May": 5, "Jun": 6,
                        "Jul": 7, "Aug": 8, "Sep": 9, "Oct": 10, "Nov": 11, "Dec": 12
                    }

                    if month_str in month_map:
                        completed = month_map[month_str]  # 月份 = 完成数
                        total = day  # 日期 = 总数
                        print(f"  [提示] Excel 日期格式 '{progress_str}' 解析为 {completed}/{total}")
                    else:
                        raise ValueError(f"无法识别月份: {month_str}")
                else:
                    raise ValueError(f"无法解析: {progress_str}")
            except Exception as e:
                print(f"  [警告] 无法解析进度字段: {progress_str} ({e})")
                return None, 0

        if completed == 0 or completed >= len(power_points):
            return None, 0

        # 读取已完成的预测数据
        predictions = np.zeros((len(power_points), 3))
        predictions[:] = np.nan  # 初始化为 NaN

        n_loaded = 0
        for i, power in enumerate(power_points[:completed]):
            w_col = f"width_um_{int(power)}W"
            d_col = f"depth_um_{int(power)}W"
            a_col = f"area_um2_{int(power)}W"

            if w_col in row and d_col in row and a_col in row:
                predictions[i, 0] = row[w_col]
                predictions[i, 1] = row[d_col]
                predictions[i, 2] = row[a_col]
                n_loaded += 1

        if n_loaded != completed:
            print(f"  [警告] 期望加载 {completed} 个功率点，实际加载 {n_loaded} 个")

        return predictions, completed

    def _load_history_csv(
        self,
        csv_path: str,
    ) -> Tuple[List[List[float]], List[float], List[Optional[np.ndarray]]]:
        """
        从 CSV 文件加载历史数据

        Returns
        -------
        tuple
            (all_x, all_y, all_predictions)
        """
        df = pd.read_csv(csv_path)

        all_x = []
        all_y = []
        all_predictions = []

        csv_col_map = get_param_csv_column_map()

        # 检查必要的列是否存在
        required_csv_cols = []
        for param_name in self.runner.active_names:
            if param_name not in csv_col_map:
                raise ValueError(f"未知参数名: {param_name}")
            csv_col = self._resolve_param_csv_column(df, param_name)
            required_csv_cols.append(csv_col)

        if not all(col in df.columns for col in required_csv_cols):
            raise ValueError(f"CSV 缺少必要的参数列: {required_csv_cols}")

        # 加载实验数据
        from ..data import load_experiment_data

        exp_data = load_experiment_data(self.config.exp_csv)
        power_points = exp_data[:, 0]
        exp_observations = exp_data[:, 1:4]
        output_weights = np.asarray(self.config.output_weights, dtype=float)

        # 检测功率点列
        power_cols_map = {}
        for power in power_points:
            w_col = f"width_um_{int(power)}W"
            d_col = f"depth_um_{int(power)}W"
            a_col = f"area_um2_{int(power)}W"
            if all(c in df.columns for c in [w_col, d_col, a_col]):
                power_cols_map[power] = (w_col, d_col, a_col)

        has_detailed_data = len(power_cols_map) == len(power_points)

        for _, row in df.iterrows():
            params = []
            for param_name in self.runner.active_names:
                csv_col = self._resolve_param_csv_column(df, param_name)
                params.append(float(row[csv_col]))

            predictions = None
            nrmse = np.nan

            # 优先使用 CSV 中已存储的 objective（热启动时保持历史目标值不被重算覆盖）
            if "objective" in df.columns and not pd.isna(row["objective"]):
                nrmse = float(row["objective"])

            # 尝试解析详细数据
            if has_detailed_data:
                try:
                    current_preds = []
                    any_nan = False
                    for power in power_points:
                        w_col, d_col, a_col = power_cols_map[power]
                        vals = [row[w_col], row[d_col], row[a_col]]
                        if any(pd.isna(v) for v in vals):
                            any_nan = True
                            break
                        current_preds.append(vals)

                    if not any_nan:
                        predictions = np.array(current_preds)
                        if np.isnan(nrmse):
                            nrmse, _ = compute_bayes_nrmse_percent(
                                predictions,
                                exp_observations,
                                weights=output_weights,
                            )
                except Exception:
                    predictions = None

            # 若 objective 缺失且无法重算，则尝试回退到旧列名 NRMSE
            if np.isnan(nrmse):
                if "NRMSE" in df.columns and not pd.isna(row["NRMSE"]):
                    nrmse = float(row["NRMSE"])

            all_x.append(params)
            all_y.append(nrmse)
            all_predictions.append(predictions)

        print(f"  从 {csv_path} 加载 {len(all_x)} 条记录")
        return all_x, all_y, all_predictions
