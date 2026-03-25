#!/usr/bin/env python3
"""
Live monitor for OpenFOAM-style solver logs.

It reparses the current log file on a timer and updates a small dashboard
covering timestep size, solver residuals, continuity errors, and melt
corrector residuals.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

TIME_RE = re.compile(r"^Time = ([^\s]+)")
DELTAT_RE = re.compile(r"^deltaT = ([^\s]+)")
SOLVER_RE = re.compile(
    r"^[A-Za-z0-9_]+:\s+Solving for ([^,]+), Initial residual = ([^,]+),"
)
CONTINUITY_RE = re.compile(
    r"^time step continuity errors : sum local = ([^,]+), "
    r"global = ([^,]+), cumulative = ([^\s]+)"
)
EPSILON_RE = re.compile(
    r"^Correcting epsilon1, mean residual = ([^,]+), max residual = ([^\s]+)"
)


def _to_float(token: str) -> float | None:
    try:
        return float(token)
    except ValueError:
        return None


def parse_log(log_path: Path) -> dict[str, object]:
    time_values: list[float] = []
    delta_t_values: list[float] = []
    continuity_step: list[int] = []
    continuity_local: list[float] = []
    continuity_global: list[float] = []
    continuity_cumulative: list[float] = []
    epsilon_step: list[int] = []
    epsilon_mean: list[float] = []
    epsilon_max: list[float] = []
    residuals: dict[str, list[tuple[int, float]]] = defaultdict(list)

    if not log_path.exists():
        return {
            "times": time_values,
            "delta_t": delta_t_values,
            "continuity_step": continuity_step,
            "continuity_local": continuity_local,
            "continuity_global": continuity_global,
            "continuity_cumulative": continuity_cumulative,
            "epsilon_step": epsilon_step,
            "epsilon_mean": epsilon_mean,
            "epsilon_max": epsilon_max,
            "residuals": residuals,
            "solver_steps": 0,
        }

    current_step = -1

    with log_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            time_match = TIME_RE.match(line)
            if time_match:
                value = _to_float(time_match.group(1))
                if value is not None:
                    time_values.append(value)
                    current_step += 1
                continue

            delta_t_match = DELTAT_RE.match(line)
            if delta_t_match:
                value = _to_float(delta_t_match.group(1))
                if value is not None:
                    delta_t_values.append(value)
                continue

            solver_match = SOLVER_RE.match(line)
            if solver_match and current_step >= 0:
                field_name = solver_match.group(1).strip()
                residual = _to_float(solver_match.group(2))
                if residual is not None and residual > 0.0:
                    residuals[field_name].append((current_step, residual))
                continue

            continuity_match = CONTINUITY_RE.match(line)
            if continuity_match and current_step >= 0:
                local_val = _to_float(continuity_match.group(1))
                global_val = _to_float(continuity_match.group(2))
                cumulative_val = _to_float(continuity_match.group(3))
                if None not in (local_val, global_val, cumulative_val):
                    continuity_step.append(current_step)
                    continuity_local.append(local_val)  # type: ignore[arg-type]
                    continuity_global.append(global_val)  # type: ignore[arg-type]
                    continuity_cumulative.append(cumulative_val)  # type: ignore[arg-type]
                continue

            epsilon_match = EPSILON_RE.match(line)
            if epsilon_match and current_step >= 0:
                mean_val = _to_float(epsilon_match.group(1))
                max_val = _to_float(epsilon_match.group(2))
                if None not in (mean_val, max_val):
                    epsilon_step.append(current_step)
                    epsilon_mean.append(mean_val)  # type: ignore[arg-type]
                    epsilon_max.append(max_val)  # type: ignore[arg-type]

    return {
        "times": time_values,
        "delta_t": delta_t_values,
        "continuity_step": continuity_step,
        "continuity_local": continuity_local,
        "continuity_global": continuity_global,
        "continuity_cumulative": continuity_cumulative,
        "epsilon_step": epsilon_step,
        "epsilon_mean": epsilon_mean,
        "epsilon_max": epsilon_max,
        "residuals": residuals,
        "solver_steps": max(current_step + 1, 0),
    }


def tail_list(values: list[float], limit: int) -> list[float]:
    if limit <= 0 or len(values) <= limit:
        return values
    return values[-limit:]


def tail_pairs(values: list[tuple[int, float]], limit: int) -> list[tuple[int, float]]:
    if limit <= 0 or len(values) <= limit:
        return values
    return values[-limit:]


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logfile", type=Path, help="Path to the solver log file")
    parser.add_argument(
        "--refresh",
        type=float,
        default=2.0,
        help="Refresh period in seconds (default: 2.0)",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=300,
        help="Maximum number of recent samples to display per curve (default: 300)",
    )
    return parser


def main() -> int:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation
    except ModuleNotFoundError as exc:
        raise SystemExit(
            "matplotlib is required for live plotting. "
            "Please run this script from your conda environment."
        ) from exc

    args = build_argument_parser().parse_args()
    log_path = args.logfile.resolve()

    fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    ax_dt, ax_res, ax_cont, ax_eps = axes.flat

    def update(_frame: int) -> None:
        data = parse_log(log_path)
        times = data["times"]
        delta_t = data["delta_t"]
        residuals = data["residuals"]

        ax_dt.clear()
        ax_res.clear()
        ax_cont.clear()
        ax_eps.clear()

        fig.suptitle(
            f"Live solver monitor: {log_path.name}",
            fontsize=14,
        )

        if times:
            shown_times = tail_list(times, args.window)
            shown_dt = tail_list(delta_t, len(shown_times))
            ax_dt.plot(shown_times[: len(shown_dt)], shown_dt, marker="o", ms=3)
            ax_dt.set_xlabel("Physical time [s]")
            ax_dt.set_ylabel("deltaT [s]")
            ax_dt.set_title("Timestep")
            ax_dt.grid(True, alpha=0.3)
        else:
            ax_dt.set_title("Timestep")
            ax_dt.text(0.5, 0.5, "Waiting for time steps...", ha="center", va="center")

        for field_name, series in sorted(residuals.items()):
            shown = tail_pairs(series, args.window)
            if shown:
                x_vals = [item[0] for item in shown]
                y_vals = [item[1] for item in shown]
                ax_res.semilogy(x_vals, y_vals, label=field_name)
        ax_res.set_title("Initial Residuals")
        ax_res.set_xlabel("Time-step index")
        ax_res.set_ylabel("Residual")
        ax_res.grid(True, which="both", alpha=0.3)
        if residuals:
            ax_res.legend(fontsize=8, ncol=2)

        cont_step = tail_list(data["continuity_step"], args.window)
        cont_local = tail_list(data["continuity_local"], args.window)
        cont_global = tail_list(data["continuity_global"], args.window)
        cont_cumulative = tail_list(data["continuity_cumulative"], args.window)
        if cont_step:
            ax_cont.semilogy(cont_step, [abs(v) for v in cont_local], label="sum local")
            ax_cont.semilogy(cont_step, [abs(v) for v in cont_global], label="global")
            ax_cont.semilogy(
                cont_step,
                [max(abs(v), sys.float_info.min) for v in cont_cumulative],
                label="cumulative",
            )
        ax_cont.set_title("Continuity Errors")
        ax_cont.set_xlabel("Time-step index")
        ax_cont.set_ylabel("|error|")
        ax_cont.grid(True, which="both", alpha=0.3)
        if cont_step:
            ax_cont.legend(fontsize=8)

        eps_step = tail_list(data["epsilon_step"], args.window)
        eps_mean = tail_list(data["epsilon_mean"], args.window)
        eps_max = tail_list(data["epsilon_max"], args.window)
        if eps_step:
            ax_eps.semilogy(eps_step, [max(v, sys.float_info.min) for v in eps_mean], label="mean")
            ax_eps.semilogy(eps_step, [max(v, sys.float_info.min) for v in eps_max], label="max")
        ax_eps.set_title("epsilon1 Corrector Residuals")
        ax_eps.set_xlabel("Time-step index")
        ax_eps.set_ylabel("Residual")
        ax_eps.grid(True, which="both", alpha=0.3)
        if eps_step:
            ax_eps.legend(fontsize=8)

        solver_steps = data["solver_steps"]
        latest_time = times[-1] if times else math.nan
        fig.texts.clear()
        fig.text(
            0.01,
            0.01,
            f"steps={solver_steps}   latestTime={latest_time if times else 'n/a'}   refresh={args.refresh:.1f}s",
            fontsize=9,
        )

    FuncAnimation(fig, update, interval=max(args.refresh, 0.2) * 1000, cache_frame_data=False)
    update(0)
    plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
