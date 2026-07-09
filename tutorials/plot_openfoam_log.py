#!/usr/bin/env python3
"""Plot residuals and run metrics from OpenFOAM-style log files.

Examples:
    python3 tutorials/plot_openfoam_log.py log.laserbeamFoam
    python3 tutorials/plot_openfoam_log.py log.laserbeamFoam --fields T p_rgh --value final
    python3 tutorials/plot_openfoam_log.py log.laserbeamFoam --metrics execution deltaT courantMax
    python3 tutorials/plot_openfoam_log.py log.laserbeamFoam --watch --refresh-interval 2 --show
    python3 tutorials/plot_openfoam_log.py log.laserbeamFoam --list-fields
"""

from __future__ import annotations

import argparse
import math
import re
import sys
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path


FLOAT_RE = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"

TIME_RE = re.compile(rf"^Time = (?P<time>{FLOAT_RE})$")
DELTA_T_RE = re.compile(rf"^deltaT = (?P<delta_t>{FLOAT_RE})$")
EXECUTION_RE = re.compile(
    rf"^ExecutionTime = (?P<execution>{FLOAT_RE}) s\s+ClockTime = (?P<clock>{FLOAT_RE}) s$"
)
COURANT_RE = re.compile(
    rf"^Courant Number mean: (?P<mean>{FLOAT_RE}) max: (?P<max>{FLOAT_RE})$"
)
INTERFACE_COURANT_RE = re.compile(
    rf"^Interface Courant Number mean: (?P<mean>{FLOAT_RE}) max: (?P<max>{FLOAT_RE})$"
)
CONTINUITY_RE = re.compile(
    rf"^time step continuity errors : sum local = (?P<sum_local>{FLOAT_RE}), "
    rf"global = (?P<global>{FLOAT_RE}), cumulative = (?P<cumulative>{FLOAT_RE})$"
)
EPSILON_RE = re.compile(
    rf"^Correcting epsilon1, mean residual = (?P<mean>{FLOAT_RE}), max residual = (?P<max>{FLOAT_RE})$"
)
SOLVE_RE = re.compile(
    rf"^(?P<solver>[^:]+):\s+Solving for (?P<field>[^,]+), "
    rf"Initial residual = (?P<initial>{FLOAT_RE}), "
    rf"Final residual = (?P<final>{FLOAT_RE}), "
    rf"No Iterations (?P<iterations>\d+)$"
)
PROGRAMMED_POSITION_RE = re.compile(
    rf"^programmed position = \((?P<x>{FLOAT_RE}) (?P<y>{FLOAT_RE}) (?P<z>{FLOAT_RE})\)$"
)
LASER_POWER_RE = re.compile(rf"^power = (?P<power>{FLOAT_RE})$")
MIN_RAY_POWER_RE = re.compile(rf"^minRayPower = (?P<min_ray_power>{FLOAT_RE})$")
RETAINED_RAYS_RE = re.compile(rf"^retained rays after gating = (?P<retained_rays>\d+)$")
SAMPLED_RAYS_RE = re.compile(
    rf"^Sampled rays for profile: (?P<sampled_rays>\d+), "
    rf"cumulative local particles: (?P<cumulative_local_particles>\d+)$"
)
TOTAL_Q_DEPOSITED_RE = re.compile(rf"^Total Q deposited: (?P<total_q>{FLOAT_RE})$")


@dataclass
class SolveEntry:
    solver: str
    field: str
    initial: float
    final: float
    iterations: int


@dataclass
class StepData:
    index: int
    time: float
    delta_t: float | None = None
    execution_time: float | None = None
    clock_time: float | None = None
    courant_mean: float | None = None
    courant_max: float | None = None
    interface_courant_mean: float | None = None
    interface_courant_max: float | None = None
    continuity_sum_local: list[float] = field(default_factory=list)
    continuity_global: list[float] = field(default_factory=list)
    continuity_cumulative: list[float] = field(default_factory=list)
    epsilon_mean: list[float] = field(default_factory=list)
    epsilon_max: list[float] = field(default_factory=list)
    solves: list[SolveEntry] = field(default_factory=list)
    laser_power: float | None = None
    min_ray_power: float | None = None
    retained_rays: int | None = None
    sampled_rays: int | None = None
    cumulative_local_particles: int | None = None
    total_q_deposited: float | None = None
    programmed_position_x: float | None = None
    programmed_position_y: float | None = None
    programmed_position_z: float | None = None


METRICS = {
    "execution": ("execution_time", "Execution time [s]", "linear"),
    "clock": ("clock_time", "Clock time [s]", "linear"),
    "deltaT": ("delta_t", "deltaT [s]", "linear"),
    "courantMean": ("courant_mean", "Courant mean [-]", "linear"),
    "courantMax": ("courant_max", "Courant max [-]", "linear"),
    "interfaceCourantMean": (
        "interface_courant_mean",
        "Interface Courant mean [-]",
        "linear",
    ),
    "interfaceCourantMax": (
        "interface_courant_max",
        "Interface Courant max [-]",
        "linear",
    ),
    "continuityLocal": ("continuity_sum_local", "Continuity sum local", "linear"),
    "continuityGlobal": ("continuity_global", "Continuity global", "linear"),
    "continuityCumulative": ("continuity_cumulative", "Continuity cumulative", "linear"),
    "epsilonMean": ("epsilon_mean", "epsilon1 mean residual", "log"),
    "epsilonMax": ("epsilon_max", "epsilon1 max residual", "log"),
    "laserPower": ("laser_power", "Laser power", "linear"),
    "minRayPower": ("min_ray_power", "Minimum ray power", "linear"),
    "retainedRays": ("retained_rays", "Retained rays after gating", "linear"),
    "sampledRays": ("sampled_rays", "Sampled rays for profile", "linear"),
    "cumulativeLocalParticles": (
        "cumulative_local_particles",
        "Cumulative local particles",
        "linear",
    ),
    "depositedEnergy": ("total_q_deposited", "Total Q deposited", "linear"),
    "programmedPositionX": ("programmed_position_x", "Programmed position x [m]", "linear"),
    "programmedPositionY": ("programmed_position_y", "Programmed position y [m]", "linear"),
    "programmedPositionZ": ("programmed_position_z", "Programmed position z [m]", "linear"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logfile", type=Path, help="Path to the log file")
    parser.add_argument(
        "--fields",
        nargs="+",
        default=None,
        help="Field names to plot, for example T p_rgh alpha.metal",
    )
    parser.add_argument(
        "--field-regex",
        default=None,
        help="Regex used to select fields from the log",
    )
    parser.add_argument(
        "--top-fields",
        type=int,
        default=5,
        help="How many most-frequent fields to plot when --fields is not given",
    )
    parser.add_argument(
        "--value",
        choices=("initial", "final", "iterations"),
        default="initial",
        help="Solver value to plot for selected fields",
    )
    parser.add_argument(
        "--reduce",
        choices=("last", "first", "max", "min", "mean"),
        default="last",
        help="How to reduce repeated solves of the same field within one time step",
    )
    parser.add_argument(
        "--x",
        choices=("step", "time", "execution"),
        default="time",
        help="X axis for the plots",
    )
    parser.add_argument(
        "--metrics",
        nargs="*",
        choices=tuple(METRICS),
        default=[],
        help="Additional non-solver metrics to plot in separate subplots",
    )
    parser.add_argument(
        "--no-residuals",
        action="store_true",
        help="Skip solver field plots and plot only metrics from --metrics",
    )
    parser.add_argument(
        "--start-time",
        type=float,
        default=None,
        help="Ignore entries before this simulation time",
    )
    parser.add_argument(
        "--end-time",
        type=float,
        default=None,
        help="Ignore entries after this simulation time",
    )
    parser.add_argument(
        "--yscale",
        choices=("auto", "linear", "log"),
        default="auto",
        help="Y scale for the residual subplot",
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Optional figure title",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Save figure to this file instead of the default <log>.png",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Open an interactive window in addition to saving",
    )
    parser.add_argument(
        "--watch",
        action="store_true",
        help="Refresh the plot while the log file is being updated",
    )
    parser.add_argument(
        "--refresh-interval",
        type=float,
        default=2.0,
        help="Seconds between refreshes when --watch is used",
    )
    parser.add_argument(
        "--list-fields",
        action="store_true",
        help="Print discovered field names and exit",
    )
    return parser.parse_args()


def parse_log(path: Path) -> list[StepData]:
    steps: list[StepData] = []
    current_step: StepData | None = None
    pending_next: dict[str, float] = {}

    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue

            match = COURANT_RE.match(line)
            if match:
                pending_next["courant_mean"] = float(match.group("mean"))
                pending_next["courant_max"] = float(match.group("max"))
                continue

            match = INTERFACE_COURANT_RE.match(line)
            if match:
                pending_next["interface_courant_mean"] = float(match.group("mean"))
                pending_next["interface_courant_max"] = float(match.group("max"))
                continue

            match = DELTA_T_RE.match(line)
            if match:
                pending_next["delta_t"] = float(match.group("delta_t"))
                continue

            match = TIME_RE.match(line)
            if match:
                current_step = StepData(index=len(steps) + 1, time=float(match.group("time")))
                for key, value in pending_next.items():
                    setattr(current_step, key, value)
                pending_next.clear()
                steps.append(current_step)
                continue

            if current_step is None:
                continue

            match = SOLVE_RE.match(line)
            if match:
                current_step.solves.append(
                    SolveEntry(
                        solver=match.group("solver").strip(),
                        field=match.group("field").strip(),
                        initial=float(match.group("initial")),
                        final=float(match.group("final")),
                        iterations=int(match.group("iterations")),
                    )
                )
                continue

            match = EXECUTION_RE.match(line)
            if match:
                current_step.execution_time = float(match.group("execution"))
                current_step.clock_time = float(match.group("clock"))
                continue

            match = PROGRAMMED_POSITION_RE.match(line)
            if match:
                current_step.programmed_position_x = float(match.group("x"))
                current_step.programmed_position_y = float(match.group("y"))
                current_step.programmed_position_z = float(match.group("z"))
                continue

            match = LASER_POWER_RE.match(line)
            if match:
                current_step.laser_power = float(match.group("power"))
                continue

            match = MIN_RAY_POWER_RE.match(line)
            if match:
                current_step.min_ray_power = float(match.group("min_ray_power"))
                continue

            match = RETAINED_RAYS_RE.match(line)
            if match:
                current_step.retained_rays = int(match.group("retained_rays"))
                continue

            match = SAMPLED_RAYS_RE.match(line)
            if match:
                current_step.sampled_rays = int(match.group("sampled_rays"))
                current_step.cumulative_local_particles = int(match.group("cumulative_local_particles"))
                continue

            match = TOTAL_Q_DEPOSITED_RE.match(line)
            if match:
                current_step.total_q_deposited = float(match.group("total_q"))
                continue

            match = CONTINUITY_RE.match(line)
            if match:
                current_step.continuity_sum_local.append(float(match.group("sum_local")))
                current_step.continuity_global.append(float(match.group("global")))
                current_step.continuity_cumulative.append(float(match.group("cumulative")))
                continue

            match = EPSILON_RE.match(line)
            if match:
                current_step.epsilon_mean.append(float(match.group("mean")))
                current_step.epsilon_max.append(float(match.group("max")))

    return steps


def reduce_values(values: list[float], mode: str) -> float | None:
    if not values:
        return None
    if mode == "last":
        return values[-1]
    if mode == "first":
        return values[0]
    if mode == "max":
        return max(values)
    if mode == "min":
        return min(values)
    if mode == "mean":
        return sum(values) / len(values)
    raise ValueError(f"unsupported reduction mode: {mode}")


def step_x(step: StepData, mode: str) -> float | None:
    if mode == "step":
        return float(step.index)
    if mode == "time":
        return step.time
    if mode == "execution":
        return step.execution_time
    raise ValueError(f"unsupported x mode: {mode}")


def filtered_steps(
    steps: list[StepData],
    start_time: float | None,
    end_time: float | None,
) -> list[StepData]:
    result: list[StepData] = []
    for step in steps:
        if start_time is not None and step.time < start_time:
            continue
        if end_time is not None and step.time > end_time:
            continue
        result.append(step)
    return result


def discovered_fields(steps: list[StepData]) -> list[str]:
    counts = Counter(solve.field for step in steps for solve in step.solves)
    return [name for name, _ in counts.most_common()]


def select_fields(args: argparse.Namespace, steps: list[StepData]) -> list[str]:
    counts = Counter(solve.field for step in steps for solve in step.solves)
    names = [name for name, _ in counts.most_common()]

    if args.fields:
        selected = [name for name in args.fields if name in counts]
        missing = [name for name in args.fields if name not in counts]
        if missing:
            print(f"warning: fields not found in log: {', '.join(missing)}", file=sys.stderr)
        return selected

    if args.field_regex:
        pattern = re.compile(args.field_regex)
        return [name for name in names if pattern.search(name)]

    return names[: max(args.top_fields, 0)]


def residual_series(
    steps: list[StepData],
    field_name: str,
    value_name: str,
    reduce_mode: str,
    x_mode: str,
) -> tuple[list[float], list[float]]:
    xs: list[float] = []
    ys: list[float] = []
    for step in steps:
        values = [getattr(solve, value_name) for solve in step.solves if solve.field == field_name]
        value = reduce_values([float(v) for v in values], reduce_mode)
        x_value = step_x(step, x_mode)
        if value is None or x_value is None:
            continue
        xs.append(x_value)
        ys.append(value)
    return xs, ys


def metric_series(
    steps: list[StepData],
    metric_key: str,
    reduce_mode: str,
    x_mode: str,
) -> tuple[list[float], list[float]]:
    attr_name = METRICS[metric_key][0]
    xs: list[float] = []
    ys: list[float] = []
    for step in steps:
        raw_value = getattr(step, attr_name)
        if isinstance(raw_value, list):
            value = reduce_values(raw_value, reduce_mode)
        else:
            value = raw_value
        x_value = step_x(step, x_mode)
        if value is None or x_value is None:
            continue
        xs.append(x_value)
        ys.append(float(value))
    return xs, ys


def default_output_path(logfile: Path, value_name: str) -> Path:
    suffix = logfile.suffix if logfile.suffix else ".log"
    return logfile.with_name(f"{logfile.name.removesuffix(suffix)}_{value_name}.png")


def choose_yscale(args: argparse.Namespace, default_scale: str) -> str:
    if args.yscale != "auto":
        return args.yscale
    return default_scale


def configure_axis(axis, x_mode: str, y_label: str, y_scale: str, title: str) -> None:
    axis.set_title(title)
    axis.set_ylabel(y_label)
    axis.set_xlabel({"step": "Time step", "time": "Simulation time [s]", "execution": "Execution time [s]"}[x_mode])
    axis.grid(True, which="both", alpha=0.3)
    if y_scale == "log":
        axis.set_yscale("log")


def sanitize_series(xs: list[float], ys: list[float], y_scale: str) -> tuple[list[float], list[float]]:
    if y_scale != "log":
        return xs, ys

    filtered_xs: list[float] = []
    filtered_ys: list[float] = []
    for x_value, y_value in zip(xs, ys):
        if y_value > 0.0 and math.isfinite(y_value):
            filtered_xs.append(x_value)
            filtered_ys.append(y_value)
    return filtered_xs, filtered_ys


def require_matplotlib():
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit(
            "matplotlib is required to plot the log. Install it with your Python environment first."
        ) from exc
    return plt


def prepare_steps(args: argparse.Namespace) -> list[StepData]:
    return filtered_steps(parse_log(args.logfile), args.start_time, args.end_time)


def render_plot(
    args: argparse.Namespace,
    steps: list[StepData],
    plt,
    figure=None,
    *,
    verbose: bool,
):
    fields = [] if args.no_residuals else select_fields(args, steps)

    if not fields and not args.metrics:
        raise ValueError("Nothing selected to plot. Use --fields, --field-regex, or --metrics.")

    subplot_count = (1 if fields else 0) + len(args.metrics)

    if figure is None:
        figure, axes = plt.subplots(subplot_count, 1, figsize=(11, 3.8 * subplot_count), squeeze=False)
    else:
        figure.clf()
        axes = figure.subplots(subplot_count, 1, squeeze=False)

    axes_flat = [row[0] for row in axes]
    axis_index = 0

    if fields:
        axis = axes_flat[axis_index]
        axis_index += 1
        residual_scale = choose_yscale(args, "linear" if args.value == "iterations" else "log")

        for field_name in fields:
            xs, ys = residual_series(steps, field_name, args.value, args.reduce, args.x)
            xs, ys = sanitize_series(xs, ys, residual_scale)
            if not xs:
                continue
            axis.plot(xs, ys, linewidth=1.5, label=field_name)

        configure_axis(
            axis,
            args.x,
            {
                "initial": "Initial residual [-]",
                "final": "Final residual [-]",
                "iterations": "Solver iterations [-]",
            }[args.value],
            residual_scale,
            f"{args.value.capitalize()} solver data ({args.reduce} per time step)",
        )
        if axis.lines:
            axis.legend(loc="best", fontsize="small")

    for metric_key in args.metrics:
        axis = axes_flat[axis_index]
        axis_index += 1
        label = METRICS[metric_key][1]
        scale = METRICS[metric_key][2]
        xs, ys = metric_series(steps, metric_key, args.reduce, args.x)
        xs, ys = sanitize_series(xs, ys, scale)
        if xs:
            axis.plot(xs, ys, linewidth=1.5)
        configure_axis(axis, args.x, label, scale, label)

    if args.title:
        figure.suptitle(args.title)

    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.98 if args.title else 1.0))

    output_path = args.output or default_output_path(args.logfile, args.value)
    figure.savefig(output_path, dpi=160, bbox_inches="tight")
    if verbose:
        print(f"Saved plot to {output_path}")

    if args.show:
        figure.canvas.draw_idle()

    return figure, output_path


def run_once(args: argparse.Namespace) -> int:
    steps = prepare_steps(args)

    if not steps:
        print("No time-step data found in the selected range.", file=sys.stderr)
        return 1

    if args.list_fields:
        for field_name in discovered_fields(steps):
            print(field_name)
        return 0

    plt = require_matplotlib()
    _, _ = render_plot(args, steps, plt, verbose=True)

    if args.show:
        plt.show()

    return 0


def log_signature(path: Path) -> tuple[int, int] | None:
    try:
        stat = path.stat()
    except FileNotFoundError:
        return None
    return stat.st_mtime_ns, stat.st_size


def run_watch(args: argparse.Namespace) -> int:
    if args.refresh_interval <= 0.0:
        print("--refresh-interval must be positive.", file=sys.stderr)
        return 1

    plt = require_matplotlib()
    if args.show:
        plt.ion()

    figure = None
    last_signature = None
    printed_wait_message = False
    print(f"Watching {args.logfile} (refresh every {args.refresh_interval:g} s)")

    try:
        while True:
            signature = log_signature(args.logfile)

            if figure is not None and args.show and not plt.fignum_exists(figure.number):
                return 0

            if signature != last_signature or figure is None:
                steps = prepare_steps(args)

                if steps:
                    if args.list_fields:
                        for field_name in discovered_fields(steps):
                            print(field_name)
                        return 0

                    figure, _ = render_plot(
                        args,
                        steps,
                        plt,
                        figure=figure,
                        verbose=last_signature is None,
                    )
                    last_signature = signature
                    printed_wait_message = False
                elif not printed_wait_message:
                    print("No time-step data found yet; waiting for log updates.", file=sys.stderr)
                    printed_wait_message = True

            if args.show:
                plt.pause(args.refresh_interval)
            else:
                time.sleep(args.refresh_interval)
    except KeyboardInterrupt:
        return 0


def main() -> int:
    args = parse_args()
    if args.watch:
        return run_watch(args)
    return run_once(args)


if __name__ == "__main__":
    raise SystemExit(main())
