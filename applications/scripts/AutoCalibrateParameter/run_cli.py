#!/usr/bin/env python3
"""
熔池参数校准主入口

统一命令行接口，支持两种优化方法：
- bayes: 贝叶斯优化 (scikit-optimize)
- gradient: 梯度优化 (scipy.optimize.least_squares)

使用示例：
    # 贝叶斯优化
    python main.py --method bayes --n-batches 50 --batch-size 5

    # 梯度优化
    python main.py --method gradient --n-directions 8 --max-nfev 300
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Optional

import numpy as np


def parse_args() -> argparse.Namespace:
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="熔池参数校准工具",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # 通用参数
    parser.add_argument(
        "--method",
        type=str,
        choices=["bayes", "gradient"],
        default="bayes",
        help="优化方法: bayes (贝叶斯优化), gradient (梯度优化) (默认: bayes)",
    )
    parser.add_argument(
        "--config",
        type=str,
        default="config.yaml",
        help="配置文件路径 (默认: config.yaml)",
    )
    parser.add_argument(
        "--exp-csv",
        type=str,
        default=None,
        help="实验数据 CSV 文件路径 (覆盖配置文件)",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="输出目录 (覆盖配置文件)",
    )
    parser.add_argument(
        "--n-proc",
        type=int,
        default=None,
        help="OpenFOAM 并行进程数 (覆盖配置文件)",
    )
    parser.add_argument(
        "--hpc",
        action="store_true",
        help="启用 HPC 模式",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="详细输出",
    )

    # 贝叶斯优化特定参数
    bayes_group = parser.add_argument_group("贝叶斯优化参数")
    bayes_group.add_argument(
        "--n-initial-points",
        type=int,
        default=None,
        help="初始采样点数 (默认: 15)",
    )
    bayes_group.add_argument(
        "--n-batches",
        type=int,
        default=None,
        help="优化批次数 (默认: 50)",
    )
    bayes_group.add_argument(
        "--batch-size",
        type=int,
        default=None,
        help="每批评估数 (默认: 5)",
    )
    bayes_group.add_argument(
        "--acq-func",
        type=str,
        choices=["EI", "PI", "LCB", "gp_hedge"],
        default=None,
        help="采集函数 (默认: EI)",
    )
    bayes_group.add_argument(
        "--warmstart-csv",
        type=str,
        default=None,
        help="历史数据 CSV 文件路径，用于热启动 (如 bayes_history.csv)",
    )

    # 梯度优化特定参数
    gradient_group = parser.add_argument_group("梯度优化参数")
    gradient_group.add_argument(
        "--n-directions",
        type=int,
        default=None,
        help="多起点方向数 (默认: 8)",
    )
    gradient_group.add_argument(
        "--max-nfev",
        type=int,
        default=None,
        help="最大函数评估次数 (默认: 300)",
    )
    gradient_group.add_argument(
        "--ftol",
        type=float,
        default=None,
        help="目标函数收敛容差 (默认: 1e-4)",
    )
    gradient_group.add_argument(
        "--xtol",
        type=float,
        default=None,
        help="参数收敛容差 (默认: 1e-4)",
    )
    gradient_group.add_argument(
        "--init-strategy",
        type=str,
        choices=["legacy_single", "legacy_multi", "directional"],
        default=None,
        help="初始化策略 (默认: legacy_multi)",
    )

    return parser.parse_args()


def create_config(args: argparse.Namespace):
    """根据命令行参数创建配置对象"""
    from AutoCalibrateParameter.config import BayesConfig, GradientConfig

    # 基础路径
    script_dir = Path(__file__).parent
    base_dir = script_dir
    
    # Check if config exists in CWD or absolute path
    if Path(args.config).exists():
        config_path = Path(args.config).resolve()
    else:
        # Fallback to script directory (but we prefer CWD for runs)
        config_path = script_dir / args.config

    # Default runs root to CWD/runs if not specified in config
    default_runs_root = Path.cwd() / "runs"

    # 根据方法选择配置类
    if args.method == "bayes":
        if config_path.exists():
            config = BayesConfig.from_yaml(config_path)
        else:
            config = BayesConfig(
                case_dir=base_dir / "openfoam_case",
                postproc_script=base_dir / "postprocess.py",
                exp_csv=base_dir / "experimental_data.csv",
                runs_root=default_runs_root,
            )

        # 覆盖命令行参数
        if args.n_initial_points is not None:
            config.n_initial_points = args.n_initial_points
        if args.n_batches is not None:
            config.n_batches = args.n_batches
        if args.batch_size is not None:
            config.batch_size = args.batch_size
        if args.acq_func is not None:
            config.acq_func = args.acq_func
        if args.warmstart_csv is not None:
            config.warmstart_csv = args.warmstart_csv

    elif args.method == "gradient":
        if config_path.exists():
            config = GradientConfig.from_yaml(config_path)
        else:
            config = GradientConfig(
                case_dir=base_dir / "openfoam_case",
                postproc_script=base_dir / "postprocess.py",
                exp_csv=base_dir / "experimental_data.csv",
                runs_root=base_dir / "runs",
            )

        # 覆盖命令行参数
        if args.n_directions is not None:
            config.n_directions = args.n_directions
        if args.max_nfev is not None:
            config.max_nfev = args.max_nfev
        if args.ftol is not None:
            config.ftol = args.ftol
        if args.xtol is not None:
            config.xtol = args.xtol
        if args.init_strategy is not None:
            config.init_strategy = args.init_strategy

    else:
        raise ValueError(f"未知方法: {args.method}")

    # 通用参数覆盖
    if args.exp_csv is not None:
        config.exp_csv = Path(args.exp_csv)
    if args.output_dir is not None:
        config.runs_root = Path(args.output_dir)
    if args.n_proc is not None:
        config.n_proc = args.n_proc
    if args.hpc:
        config.hpc_mode = True

    return config


def run_bayes(config, exp_data: np.ndarray, verbose: bool = False):
    """运行贝叶斯优化"""
    from AutoCalibrateParameter.optimizers import BayesOptimizer
    from AutoCalibrateParameter.simulation import SimulationRunner

    print("\n" + "=" * 60)
    print("贝叶斯优化 (scikit-optimize)")
    print("=" * 60)
    print(f"初始点数: {config.n_initial_points}")
    print(f"批次数: {config.n_batches}")
    print(f"每批大小: {config.batch_size}")
    print(f"采集函数: {config.acq_func}")
    print()

    runner = SimulationRunner(config)
    optimizer = BayesOptimizer(config, runner)
    result = optimizer.optimize(exp_data)

    return result


def run_gradient(config, exp_data: np.ndarray, verbose: bool = False):
    """运行梯度优化"""
    from AutoCalibrateParameter.optimizers import GradientOptimizer
    from AutoCalibrateParameter.simulation import SimulationRunner

    print("\n" + "=" * 60)
    print("梯度优化 (scipy.optimize.least_squares)")
    print("=" * 60)
    print(f"方向数: {config.n_directions}")
    print(f"最大评估: {config.max_nfev}")
    print(f"初始化策略: {config.init_strategy}")
    print()

    runner = SimulationRunner(config)
    optimizer = GradientOptimizer(config, runner)
    result = optimizer.optimize(exp_data)

    return result


def print_result(result, param_names: list = None):
    """打印优化结果"""
    if param_names is None:
        from AutoCalibrateParameter.utils import get_param_names

        param_names = get_param_names()

    print("\n" + "=" * 60)
    print("优化结果")
    print("=" * 60)
    print(f"方法: {result.method}")
    print(f"最优目标值: {result.best_cost:.6e}")
    print(f"总评估次数: {result.n_evaluations}")
    print(f"消息: {result.message}")
    print()
    print("最优参数:")
    for name, val in zip(param_names, result.best_params):
        print(f"  {name}: {val:.6g}")


def main():
    """主函数"""
    args = parse_args()

    print("=" * 60)
    print("熔池参数校准工具")
    print("=" * 60)
    print(f"方法: {args.method}")
    print(f"配置: {args.config}")
    print()

    # 创建配置
    try:
        config = create_config(args)
    except Exception as e:
        print(f"错误: 无法创建配置 - {e}")
        sys.exit(1)

    # 加载实验数据
    from AutoCalibrateParameter.data import load_experiment_data

    try:
        exp_data = load_experiment_data(config.exp_csv)
        print(f"加载实验数据: {config.exp_csv}")
        print(f"  数据点数: {len(exp_data)}")
        print(f"  功率范围: {exp_data[:, 0].min():.0f} - {exp_data[:, 0].max():.0f} W")
    except Exception as e:
        print(f"错误: 无法加载实验数据 - {e}")
        sys.exit(1)

    # 创建输出目录
    # 直接使用 runs_root，不再自动追加 "{method}_optimization" 子目录
    output_dir = config.runs_root

    # 兼容旧配置：若 runs_root 已是类似 ".../bayes_optimization"，自动回退到上一级 runs
    legacy_method_dir = f"{args.method}_optimization"
    if output_dir.name == legacy_method_dir:
        output_dir = output_dir.parent

    output_dir.mkdir(parents=True, exist_ok=True)

    # 重要：更新 config.runs_root 以确保 CSV 保存到正确位置
    config.runs_root = output_dir
    print(f"输出目录: {output_dir}")

    # 保存配置
    config_save_path = output_dir / "config.yaml"
    config.save_yaml(config_save_path)
    print(f"配置已保存: {config_save_path}")

    # 运行优化
    start_time = time.time()

    try:
        if args.method == "bayes":
            result = run_bayes(config, exp_data, args.verbose)
        elif args.method == "gradient":
            result = run_gradient(config, exp_data, args.verbose)
        else:
            raise ValueError(f"未知方法: {args.method}")
    except KeyboardInterrupt:
        print("\n用户中断")
        sys.exit(1)
    except Exception as e:
        print(f"\n错误: 优化失败 - {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    elapsed = time.time() - start_time
    print(f"\n总耗时: {elapsed / 60:.1f} 分钟")

    # 打印结果
    print_result(result)

    # 保存结果
    result_path = output_dir / "result.json"
    result.save(result_path)
    print(f"\n结果已保存: {result_path}")

    # 生成对比图
    try:
        from AutoCalibrateParameter.postprocess import (
            save_comparison_plot,
            plot_predictions_vs_experiments,
        )
        from AutoCalibrateParameter.simulation import SimulationRunner

        runner = SimulationRunner(config)
        power_points = exp_data[:, 0]
        observations = exp_data[:, 1:4]
        predictions = runner.run(result.best_params, power_points)

        plots_dir = output_dir / "plots"
        plots_dir.mkdir(exist_ok=True)

        save_comparison_plot(
            plots_dir,
            power_points,
            observations,
            predictions,
            nrmse=result.best_cost,
            title=f"{result.method} 优化结果",
        )

        plot_predictions_vs_experiments(
            plots_dir,
            observations,
            predictions,
            title=f"{result.method} 预测 vs 实验",
        )

        print(f"图表已保存: {plots_dir}")

    except Exception as e:
        print(f"警告: 无法生成图表 - {e}")

    print("\n完成!")
    return 0


if __name__ == "__main__":
    sys.exit(main())
