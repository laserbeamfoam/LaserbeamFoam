#!/usr/bin/env python3
"""Measure completed cases with SimToPC and update results/ and validation.png.

    python3 addResults.py [case ...]

Each case is a folder of cases/ that has been run (./Allrun_long <case>), so
that it holds its last time directory with meltHistory and alpha.metal.
Without arguments, every such case is measured. SimToPC must be installed
(command "simtopc"; another command can be given in the environment variable
SIMTOPC).

What it does:
  1. reads the scan speed, the power and the spot size of each case from its
     files, and measures the cases of each scan speed together with the
     metallographic mode of SimToPC;
  2. writes the row of each case to results/simtopcSummary.csv and its
     profile to results/<case>/metallographic_profile.csv;
  3. runs plotValidation.py, which redraws validation.png.

The window: for a power and scan speed that results/simtopcSummary.csv does
not have yet, SimToPC finds the window from the steady stretches of the
cases measured. For a power and scan speed that it already has, the window
of the existing cases is kept, so that their numbers do not change and all
the cases of that speed are measured over the same part of the track; a
warning is printed if that window is not inside the steady stretch of a
new case.

Only the small files of results/ are meant to be committed, not the fields.
"""
import csv
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SUMMARY = HERE / "results" / "simtopcSummary.csv"
NUMBER = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
CONFIG = """\
mesh_density: cases
parameters_file: parameters.txt
output_dir: results
running_on: local

measure:
  y_begin: {begin}
  y_end:   {end}
  x_min:   {xMin}
  x_max:   {xMax}
  cell_size: auto

  geometry_mode: metallographic
  axes: {{track: z, vertical: y}}
  substrate_top: {plateTop}
  melt_field: meltHistory
  metal_field: alpha.metal
  width_max_rows_down: 2

  track_reference: laser_start
  window_mode: {mode}
  steady_tolerance: 5e-6
  steady_filter_length: 50e-6
  window_margin: 50e-6
  window_step: 50e-6
  window_length: 400e-6
"""


def uncommented(path):
    text = re.sub(r"/\*.*?\*/", "", path.read_text(), flags=re.S)
    return re.sub(r"//.*", "", text)


def lastTime(case):
    times = [d for d in case.iterdir() if d.is_dir() and re.fullmatch(NUMBER, d.name) and float(d.name) > 0]
    if not times:
        return None
    t = max(times, key=lambda d: float(d.name))
    return t if (t / "meltHistory").exists() and (t / "alpha.metal").exists() else None


def describe(case):
    """Scan speed (m/s), power (W), spot size (m), track length (m), plate top (m) and x range of a case."""
    pos = [tuple(map(float, m)) for m in re.findall(
        rf"\(\s*({NUMBER})\s+\(\s*({NUMBER})\s+({NUMBER})\s+({NUMBER})\s*\)\s*\)",
        uncommented(case / "constant" / "timeVsLaserPosition"))]
    power = [tuple(map(float, m)) for m in re.findall(
        rf"\(\s*({NUMBER})\s+({NUMBER})\s*\)", uncommented(case / "constant" / "timeVsLaserPower"))]
    # The laser must be on from the first to the last point of its path, which must be one straight
    # segment along z: the speed and the length of the track are those of that segment
    tOn, tOff = min(t for t, p in power if p > 0), max(t for t, p in power if p > 0)
    if len(pos) != 2 or pos[0][0] != tOn or pos[1][0] != tOff or pos[0][1:3] != pos[1][1:3] or pos[0][3] == pos[1][3]:
        sys.exit(f"{case.name}: the laser path must be a single straight segment along z, from the time the "
                 f"laser is switched on ({tOn:g} s) to the time it is switched off ({tOff:g} s); "
                 "constant/timeVsLaserPosition and constant/timeVsLaserPower do not say so")
    (t0, _, _, z0), (t1, _, _, z1) = pos
    speed = abs(z1 - z0) / (t1 - t0)
    radius = float(re.search(rf"laserRadius\s+({NUMBER})", uncommented(case / "constant" / "LaserProperties")).group(1))
    plate = uncommented(case / "system" / "bedPlateDict")
    value = lambda key: float(re.search(rf"\b{key}\s+({NUMBER})", plate).group(1))
    return dict(speed=round(speed, 6), power=max(p for _, p in power), spot=2 * radius, track=speed * (tOff - tOn),
                plateTop=value("ymax"), xMin=value("xmin"), xMax=value("xmax"))


def simtopc(names, info, mode, begin, end):
    """Summary rows and profile files of the cases names, measured together."""
    command = shlex.split(os.environ.get("SIMTOPC", "simtopc"))
    if shutil.which(command[0]) is None:
        sys.exit(f"'{command[0]}' not found: install SimToPC (https://github.com/laserbeamfoam/SimToPC)")
    work = Path(tempfile.mkdtemp(prefix="simtopc_", dir=HERE))
    try:
        first = info[names[0]]
        (work / "config.yml").write_text(CONFIG.format(begin=begin, end=end, mode=mode, plateTop=first["plateTop"],
                                                       xMin=first["xMin"], xMax=first["xMax"]))
        with open(work / "parameters.txt", "w") as f:
            f.write("# Scanning Speed (m/s) | Power (W) | Spot size (m)\n")
            for n in names:
                f.write(f"  {info[n]['speed']:g}   {info[n]['power']:g}   {info[n]['spot']:g}\n")
        for i, n in enumerate(names, 1):
            case = HERE / "cases" / n
            link = work / "cases" / f"test_case_{i}"
            link.mkdir(parents=True)
            for item in (case / "system", case / "constant", lastTime(case)):
                os.symlink(item.resolve(), link / item.name)
        log = work / "log.simtopc"
        with open(log, "w") as f:
            done = subprocess.run(command + ["measure", "config.yml"], cwd=work, stdout=f, stderr=subprocess.STDOUT)
        if done.returncode:
            sys.exit(f"SimToPC failed on {', '.join(names)}:\n" + "".join(log.read_text().splitlines(True)[-15:]))
        with open(work / "metallographic_summary.csv") as f:
            rows = list(csv.DictReader(f))
        if len(rows) != len(names):
            sys.exit(f"SimToPC gave {len(rows)} rows for the {len(names)} cases {', '.join(names)}")
        profiles = []
        for i, (n, row) in enumerate(zip(names, rows), 1):
            row["case"] = n
            profiles.append((work / "cases" / f"test_case_{i}" / "measure_results"
                             / "metallographic_profile.csv").read_text())
        return rows, profiles
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main():
    completed = sorted(c.name for c in (HERE / "cases").iterdir() if c.is_dir() and lastTime(c))
    names = sys.argv[1:] or completed
    missing = [n for n in names if n not in completed]
    if missing or not names:
        sys.exit("no completed run (last time directory with meltHistory and alpha.metal) in cases/: "
                 + (", ".join(missing) or "none found"))
    info = {n: describe(HERE / "cases" / n) for n in names}

    with open(SUMMARY) as f:
        reader = csv.DictReader(f)
        columns, old = reader.fieldnames, list(reader)
    kept = [r for r in old if r["case"] not in names]
    newRows, newProfiles = {}, {}
    for key in sorted({(info[n]["power"], info[n]["speed"]) for n in names}, reverse=True):
        group = [n for n in names if (info[n]["power"], info[n]["speed"]) == key]
        for a in ("plateTop", "xMin", "xMax"):
            if len({info[n][a] for n in group}) > 1:
                sys.exit(f"the cases {', '.join(group)} do not have the same {a}")
        search = min(info[n]["track"] for n in group)
        # SimToPC itself leaves out of the choice of the window the cases with cells larger than
        # steady_tolerance (5 um): they are measured in the window that the other cases set
        rows, profiles = simtopc(group, info, "steady_stretch", 0, search)
        window = (float(rows[0]["window_begin"]), float(rows[0]["window_end"]))
        same = [r for r in kept if (float(r["power"]), float(r["scan_speed"])) == key]
        if same:
            existing = (float(same[0]["window_begin"]), float(same[0]["window_end"]))
            if max(abs(existing[0] - window[0]), abs(existing[1] - window[1])) > 1e-9:
                stretches = {r["case"]: (r["steady_stretch_begin"], r["steady_stretch_end"]) for r in rows}
                first = profiles
                rows, profiles = simtopc(group, info, "fixed", existing[0], existing[1])
                # the profile of a fixed window has no steady_stretch column: it is taken from the first run,
                # which must have the same sections
                joined = []
                for n, p, q in zip(group, profiles, first):
                    lines, before = p.splitlines(), q.splitlines()
                    if [l.split(",", 2)[:2] for l in lines] != [l.split(",", 2)[:2] for l in before]:
                        sys.exit(f"{n}: the profiles of the two SimToPC runs do not have the same sections "
                                 f"({len(before) - 1} and {len(lines) - 1}); nothing was written")
                    joined.append("".join(f"{l},{b.rsplit(',', 1)[1]}\n" for l, b in zip(lines, before)))
                profiles = joined
                for r in rows:
                    r["steady_stretch_begin"], r["steady_stretch_end"] = stretches[r["case"]]
                    if not (float(stretches[r["case"]][0]) <= existing[0] + 1e-9
                            and existing[1] <= float(stretches[r["case"]][1]) + 1e-9):
                        print(f"WARNING: the window {1e6*existing[0]:.0f}-{1e6*existing[1]:.0f} um of the existing "
                              f"cases of {key[1]:g} m/s is not inside the steady stretch of {r['case']} "
                              f"({1e6*float(stretches[r['case']][0]):.1f}-{1e6*float(stretches[r['case']][1]):.1f} um)")
            window = existing
        print(f"{key[0]:g} W, {key[1]:g} m/s: window {1e6*window[0]:.0f}-{1e6*window[1]:.0f} um from the laser start"
              + (" (that of the existing cases)" if same else ""))
        for r, p in zip(rows, profiles):
            newRows[r["case"]], newProfiles[r["case"]] = r, p

    order = [r["case"] for r in old] + [n for n in names if n not in {r["case"] for r in old}]
    byName = {r["case"]: r for r in kept}
    byName.update(newRows)
    # a column that only a newer SimToPC writes is added (empty in the rows kept); the file is replaced
    # only once it has been written completely
    columns = columns + [k for r in newRows.values() for k in r if k not in columns and k is not None]
    columns = list(dict.fromkeys(columns))
    partial = SUMMARY.with_suffix(".csv.partial")
    with open(partial, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, restval="", lineterminator="\n")
        writer.writeheader()
        writer.writerows(byName[n] for n in order)
    os.replace(partial, SUMMARY)
    for n, text in newProfiles.items():
        (HERE / "results" / n).mkdir(exist_ok=True)
        (HERE / "results" / n / "metallographic_profile.csv").write_text(text)
        r = newRows[n]
        print(f"  {n}: width {1e6*float(r['width_mean']):.1f} +- {1e6*float(r['width_sd']):.1f} um, "
              f"depth {1e6*float(r['depth_mean']):.1f} +- {1e6*float(r['depth_sd']):.1f} um, "
              f"fill {100*float(r['fill_mean']):.1f} %")
    subprocess.run([sys.executable, str(HERE / "plotValidation.py")], check=True)


if __name__ == "__main__":
    main()
