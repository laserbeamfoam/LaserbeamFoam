#!/usr/bin/env python3

from __future__ import annotations

import argparse
import concurrent.futures
import html
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path


CASES = [
    {
        "path": "laserbeamFoam/bessel_beam_example",
        "smoke_supported": True,
    },
    {
        "path": "laserbeamFoam/core_ring_example",
        "smoke_supported": True,
    },
    {
        "path": "laserbeamFoam/csv_example",
        "smoke_supported": True,
    },
    {
        "path": "laserbeamFoam/noise_modes_example",
        "smoke_supported": True,
    },
    {
        "path": "laserbeamFoam/opa_cbc_example",
        "smoke_supported": True,
    },
    {
        "path": "laserbeamFoam/supergaussian_example",
        "smoke_supported": True,
    },
    {
        "path": "laserbeamFoam/cailabs_ringCore/core_ref",
        "enabled": False,
        "skip_reason": "temporarily excluded while the cailabs_ringCore setup is still being refined",
    },
    {
        "path": "laserbeamFoam/cailabs_ringCore/ring_core",
        "enabled": False,
        "skip_reason": "temporarily excluded while the cailabs_ringCore setup is still being refined",
    },
    {
        "path": "compressiblelaserbeamFoam/LPBF_small_vapour",
        "smoke_supported": False,
        "smoke_reason": (
            "known unstable as a shortened smoke run because the copied case "
            "can reach a zero-deposition first step and trip an FPE"
        ),
    },
    {
        "path": "compressiblelaserbeamFoam/multiComponentLaserIrradiation",
        "smoke_supported": False,
        "smoke_reason": (
            "known unstable as a shortened smoke run because the compressible "
            "multiphase initialization can trigger an early FPE before the "
            "profile case reaches a representative state"
        ),
    },
]

TIME_RE = re.compile(r"^Time = ([^\\s]+)")
Q_RE = re.compile(r"Total Q deposited:\\s*([^\\s]+)")
COURANT_RE = re.compile(r"^Courant Number mean:\\s*([^\\s]+) max:\\s*([^\\s]+)")
EXEC_RE = re.compile(r"ExecutionTime =\\s*([^\\s]+) s")


def run_command(cmd: list[str], cwd: Path | None = None) -> None:
    subprocess.run(cmd, cwd=cwd, check=True)


def foam_dictionary_set(target: Path, entry: str, value: str) -> None:
    foam_dictionary = shutil.which("foamDictionary")
    if foam_dictionary is None or not target.exists():
        return

    subprocess.run(
        [foam_dictionary, "-entry", entry, "-set", value, str(target)],
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def rewrite_for_smoke(case_dir: Path, np: int) -> None:
    control_dict = case_dir / "system" / "controlDict"
    decompose_dict = case_dir / "system" / "decomposeParDict"

    foam_dictionary_set(control_dict, "startFrom", "startTime")
    foam_dictionary_set(control_dict, "startTime", "0")
    foam_dictionary_set(control_dict, "endTime", "5e-7")
    foam_dictionary_set(control_dict, "deltaT", "1e-7")
    foam_dictionary_set(control_dict, "maxDeltaT", "5e-7")
    foam_dictionary_set(control_dict, "writeInterval", "5e-7")
    foam_dictionary_set(decompose_dict, "numberOfSubdomains", str(np))


def copy_case(root: Path, tmp_root: Path, case_rel: str, mode: str, np: int) -> Path:
    src = root / case_rel
    dst = tmp_root / case_rel

    if dst.exists():
        shutil.rmtree(dst)

    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(src, dst)

    if mode == "smoke":
        rewrite_for_smoke(dst, np)

    return dst


def parse_float(text: str) -> float | None:
    try:
        return float(text)
    except ValueError:
        return None


def commit_step(
    steps: list[dict[str, float | None]],
    current_time: float | None,
    current_q: float | None,
    current_co_max: float | None,
    current_exec: float | None,
) -> None:
    if current_time is None:
        return

    steps.append(
        {
            "time": current_time,
            "q": current_q,
            "co_max": current_co_max,
            "exec": current_exec,
        }
    )


def parse_solver_log(log_path: Path) -> list[dict[str, float | None]]:
    steps: list[dict[str, float | None]] = []
    current_time: float | None = None
    current_q: float | None = None
    current_co_max: float | None = None
    current_exec: float | None = None

    try:
        with log_path.open("r", encoding="ascii", errors="ignore") as handle:
            for raw_line in handle:
                line = raw_line.strip()

                match = TIME_RE.match(line)
                if match:
                    commit_step(steps, current_time, current_q, current_co_max, current_exec)
                    current_time = parse_float(match.group(1))
                    current_q = None
                    current_co_max = None
                    current_exec = None
                    continue

                match = Q_RE.search(line)
                if match:
                    current_q = parse_float(match.group(1))
                    continue

                match = COURANT_RE.match(line)
                if match:
                    current_co_max = parse_float(match.group(2))
                    continue

                match = EXEC_RE.search(line)
                if match:
                    current_exec = parse_float(match.group(1))
                    commit_step(steps, current_time, current_q, current_co_max, current_exec)
                    current_time = None
                    current_q = None
                    current_co_max = None
                    current_exec = None
    except FileNotFoundError:
        return []

    commit_step(steps, current_time, current_q, current_co_max, current_exec)
    return steps


def solver_logs_for_case(run_dir: Path) -> list[Path]:
    return sorted(
        path
        for path in run_dir.glob("log.*")
        if path.name not in {"log.AllrunSuite", "log.blockMesh", "log.setFields", "log.setSolidFraction", "log.decomposePar", "log.transformPoints"}
    )


def best_solver_log(run_dir: Path) -> Path | None:
    best_path: Path | None = None
    best_steps = -1

    for candidate in solver_logs_for_case(run_dir):
        steps = len(parse_solver_log(candidate))
        if steps > best_steps:
            best_steps = steps
            best_path = candidate

    return best_path


def fmt_value(value: float | None, fmt: str = ".3g") -> str:
    if value is None:
        return "-"
    return format(value, fmt)


def metric_svg(
    steps: list[dict[str, float | None]],
    metric_key: str,
    color: str,
    title: str,
    width: int = 320,
    height: int = 110,
) -> str:
    values = [
        (step["time"], step[metric_key])
        for step in steps
        if step["time"] is not None and step[metric_key] is not None
    ]

    if not values:
        return (
            f"<svg viewBox='0 0 {width} {height}' width='{width}' height='{height}'>"
            f"<rect x='0' y='0' width='{width}' height='{height}' fill='#fffdf5' stroke='#d8d1c0'/>"
            f"<text x='12' y='20' font-size='12' fill='#6f6657'>{html.escape(title)}: no data yet</text>"
            "</svg>"
        )

    pad_x = 16
    pad_y = 14
    chart_w = width - 2 * pad_x
    chart_h = height - 2 * pad_y
    x_min = min(point[0] for point in values)
    x_max = max(point[0] for point in values)
    y_min = min(point[1] for point in values)
    y_max = max(point[1] for point in values)

    if x_max <= x_min:
        x_max = x_min + 1.0
    if y_max <= y_min:
        y_max = y_min + 1.0

    coords = []
    for x_val, y_val in values:
        px = pad_x + chart_w * (x_val - x_min) / (x_max - x_min)
        py = pad_y + chart_h * (1.0 - (y_val - y_min) / (y_max - y_min))
        coords.append(f"{px:.2f},{py:.2f}")

    latest = values[-1][1]
    return (
        f"<svg viewBox='0 0 {width} {height}' width='{width}' height='{height}'>"
        f"<rect x='0' y='0' width='{width}' height='{height}' fill='#fffdf5' stroke='#d8d1c0'/>"
        f"<text x='12' y='20' font-size='12' fill='#3c3428'>{html.escape(title)}</text>"
        f"<text x='{width - 12}' y='20' font-size='12' text-anchor='end' fill='#3c3428'>{html.escape(fmt_value(latest, '.4g'))}</text>"
        f"<line x1='{pad_x}' y1='{height - pad_y}' x2='{width - pad_x}' y2='{height - pad_y}' stroke='#d8d1c0'/>"
        f"<line x1='{pad_x}' y1='{pad_y}' x2='{pad_x}' y2='{height - pad_y}' stroke='#d8d1c0'/>"
        f"<polyline fill='none' stroke='{color}' stroke-width='2' points='{' '.join(coords)}'/>"
        "</svg>"
    )


def render_dashboard(
    dashboard_path: Path,
    tmp_root: Path,
    case_statuses: dict[str, str],
    case_logs: dict[str, str],
) -> None:
    cards: list[str] = []

    for case_path in sorted(case_statuses):
        run_dir = tmp_root / case_path
        solver_log = best_solver_log(run_dir)
        steps = parse_solver_log(solver_log) if solver_log else []
        latest = steps[-1] if steps else {}
        status = html.escape(case_statuses[case_path])
        suite_log = html.escape(case_logs.get(case_path, str(run_dir / "log.AllrunSuite")))
        solver_log_label = html.escape(solver_log.name if solver_log else "-")
        solver_log_path = html.escape(str(solver_log) if solver_log else "-")

        cards.append(
            (
                "<section class='card'>"
                f"<h2>{html.escape(case_path)}</h2>"
                f"<p class='meta'><strong>Status:</strong> {status}</p>"
                f"<p class='meta'><strong>Latest time:</strong> {html.escape(fmt_value(latest.get('time') if latest else None, '.4g'))} s</p>"
                f"<p class='meta'><strong>Latest Q:</strong> {html.escape(fmt_value(latest.get('q') if latest else None, '.4g'))}</p>"
                f"<p class='meta'><strong>Latest max Co:</strong> {html.escape(fmt_value(latest.get('co_max') if latest else None, '.4g'))}</p>"
                f"<p class='meta'><strong>Suite log:</strong> {suite_log}</p>"
                f"<p class='meta'><strong>Solver log:</strong> {solver_log_label}</p>"
                f"<p class='meta small'>{solver_log_path}</p>"
                f"{metric_svg(steps, 'q', '#c65d2e', 'Total Q Deposited')}"
                f"{metric_svg(steps, 'co_max', '#1f6aa5', 'Max Courant Number')}"
                "</section>"
            )
        )

    dashboard_path.write_text(
        (
            "<!doctype html><html><head><meta charset='utf-8'>"
            "<meta http-equiv='refresh' content='2'>"
            "<title>Beam Profile Suite Monitor</title>"
            "<style>"
            "body{font-family:Helvetica,Arial,sans-serif;background:#f4efe4;color:#2f2419;margin:24px;}"
            "h1{margin:0 0 8px 0;font-size:28px;} p.top{margin:0 0 20px 0;}"
            ".grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:18px;}"
            ".card{background:#fffaf0;border:1px solid #d8d1c0;padding:16px;box-shadow:0 3px 10px rgba(0,0,0,0.06);}"
            ".card h2{margin:0 0 10px 0;font-size:18px;word-break:break-word;}"
            ".meta{margin:4px 0;font-size:13px;}"
            ".small{font-size:11px;color:#6f6657;word-break:break-all;}"
            "svg{display:block;width:100%;height:auto;margin-top:10px;}"
            "</style></head><body>"
            "<h1>Beam-Profile Suite Monitor</h1>"
            f"<p class='top'>Copied cases and outputs are under {html.escape(str(tmp_root))}. This page refreshes every 2 seconds.</p>"
            "<div class='grid'>"
            + "".join(cards)
            + "</div></body></html>"
        ),
        encoding="utf-8",
    )


def monitor_dashboard(
    dashboard_path: Path,
    tmp_root: Path,
    case_statuses: dict[str, str],
    case_logs: dict[str, str],
    lock: threading.Lock,
    stop_event: threading.Event,
) -> None:
    while not stop_event.is_set():
        with lock:
            statuses = dict(case_statuses)
            logs = dict(case_logs)
        render_dashboard(dashboard_path, tmp_root, statuses, logs)
        stop_event.wait(2.0)

    with lock:
        statuses = dict(case_statuses)
        logs = dict(case_logs)
    render_dashboard(dashboard_path, tmp_root, statuses, logs)


def run_case(root: Path, tmp_root: Path, case_rel: str, mode: str, np: int) -> None:
    run_dir = copy_case(root, tmp_root, case_rel, mode, np)
    wrapper_log = run_dir / "log.AllrunSuite"
    print(
        f"==> Running {case_rel} ({mode}) in {run_dir}"
        f" [log: {wrapper_log.name}]"
    )

    env = os.environ.copy()
    env.setdefault("FOAM_RUN_NP", str(np))
    with wrapper_log.open("w", encoding="ascii", errors="ignore") as log_file:
        subprocess.run(
            [
                "bash",
                "-lc",
                (
                    "shopt -s expand_aliases\n"
                    "alias of2506='source /lib/openfoam/openfoam2506/etc/bashrc'\n"
                    "of2506 && "
                    "bash ./Allrun"
                ),
            ],
            cwd=run_dir,
            env=env,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            check=True,
        )


def selected_cases(filter_text: str) -> list[dict[str, object]]:
    return [case for case in CASES if filter_text in str(case["path"])]


def run_case_task(
    root: Path,
    tmp_root: Path,
    case_path: str,
    mode: str,
    np: int,
) -> tuple[str, bool, str]:
    log_path = str(tmp_root / case_path / "log.AllrunSuite")

    try:
        run_case(root, tmp_root, case_path, mode, np)
        return case_path, True, log_path
    except subprocess.CalledProcessError:
        return case_path, False, log_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the beam-profile tutorial suite from copied cases under /tmp "
            "so the tracked tutorial directories stay clean."
        )
    )
    parser.add_argument(
        "--mode",
        choices=("smoke", "full"),
        default="full",
        help="Run the short checked-in case settings as-is, or shorten them even further.",
    )
    parser.add_argument(
        "--filter",
        default="",
        help="Run only cases whose relative path contains this substring.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="Print the case list and exit.",
    )
    parser.add_argument(
        "--tmp-root",
        default=os.path.join(os.environ.get("TMPDIR", "/tmp"), "beam_profiles_runs"),
        help="Destination for copied run directories.",
    )
    parser.add_argument(
        "--np",
        type=int,
        default=int(os.environ.get("FOAM_RUN_NP", "4")),
        help="Parallel process count passed through to copied case runs.",
    )
    parser.add_argument(
        "--include-unstable-smoke",
        action="store_true",
        help="Include cases that are marked as unsupported in smoke mode.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of tutorial cases to run concurrently.",
    )
    parser.add_argument(
        "--no-live-plot",
        action="store_true",
        help="Disable the live HTML dashboard under the run directory.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = Path(__file__).resolve().parent
    tmp_root = Path(args.tmp_root)
    dashboard_path = tmp_root / "live_dashboard.html"

    if args.list:
        for case in CASES:
            print(case["path"])
        return 0

    tmp_root.mkdir(parents=True, exist_ok=True)
    print(f"Using copied cases under {tmp_root}")
    if not args.no_live_plot:
        print(f"Live dashboard: {dashboard_path}")

    selected = selected_cases(args.filter)
    if not selected:
        print("No matching cases found.", file=sys.stderr)
        return 1

    if args.jobs < 1:
        print("--jobs must be at least 1.", file=sys.stderr)
        return 1

    failures: list[str] = []
    skipped: list[str] = []
    runnable: list[str] = []
    case_statuses: dict[str, str] = {}
    case_logs: dict[str, str] = {}
    status_lock = threading.Lock()
    stop_event = threading.Event()

    for case in selected:
        case_path = str(case["path"])
        case_statuses[case_path] = "queued"
        case_logs[case_path] = str(tmp_root / case_path / "log.AllrunSuite")

    for case in selected:
        case_path = str(case["path"])
        enabled = bool(case.get("enabled", True))
        smoke_supported = bool(case.get("smoke_supported", True))

        if not enabled:
            reason = str(case.get("skip_reason", "excluded from the suite"))
            print(f"==> Skipping {case_path}: {reason}")
            skipped.append(case_path)
            case_statuses[case_path] = "skipped"
            continue

        if (
            args.mode == "smoke"
            and not smoke_supported
            and not args.include_unstable_smoke
        ):
            reason = str(case.get("smoke_reason", "unsupported in smoke mode"))
            print(f"==> Skipping {case_path} (smoke): {reason}")
            skipped.append(case_path)
            case_statuses[case_path] = "skipped"
            continue

        runnable.append(case_path)

    monitor_thread: threading.Thread | None = None
    if not args.no_live_plot:
        monitor_thread = threading.Thread(
            target=monitor_dashboard,
            args=(dashboard_path, tmp_root, case_statuses, case_logs, status_lock, stop_event),
            daemon=True,
        )
        monitor_thread.start()

    try:
        if not runnable:
            pass
        elif args.jobs == 1:
            for case_path in runnable:
                with status_lock:
                    case_statuses[case_path] = "running"
                case_name, ok, log_path = run_case_task(
                    root, tmp_root, case_path, args.mode, args.np
                )
                with status_lock:
                    case_statuses[case_name] = "passed" if ok else "failed"
                if not ok:
                    failures.append(case_name)
                    print(f"==> FAILED {case_name} [log: {log_path}]", file=sys.stderr)
        else:
            max_workers = min(args.jobs, len(runnable))
            print(f"==> Running up to {max_workers} tutorial cases concurrently")

            with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_to_case = {}
                for case_path in runnable:
                    with status_lock:
                        case_statuses[case_path] = "running"
                    future = executor.submit(
                        run_case_task,
                        root,
                        tmp_root,
                        case_path,
                        args.mode,
                        args.np,
                    )
                    future_to_case[future] = case_path

                for future in concurrent.futures.as_completed(future_to_case):
                    case_name, ok, log_path = future.result()
                    with status_lock:
                        case_statuses[case_name] = "passed" if ok else "failed"
                    if not ok:
                        failures.append(case_name)
                        print(f"==> FAILED {case_name} [log: {log_path}]", file=sys.stderr)
    finally:
        if monitor_thread is not None:
            stop_event.set()
            monitor_thread.join()

    print()
    print("Suite summary")
    print(f"  passed: {len(selected) - len(failures) - len(skipped)}")
    print(f"  skipped: {len(skipped)}")
    print(f"  failed: {len(failures)}")

    if skipped:
        print("  skipped cases:")
        for case_path in skipped:
            print(f"    {case_path}")

    if failures:
        print("  failed cases:", file=sys.stderr)
        for case_path in failures:
            print(f"    {case_path}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
