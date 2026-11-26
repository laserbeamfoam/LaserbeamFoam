#!/usr/bin/env python3
"""Melt-pool width and depth of the cases, measured with SimToPC
(results/simtopcSummary.csv), against the experiments of Dilip et al. (2017)
(results/dilip2017.csv). Prints the comparison and writes validation.png.

    python3 plotValidation.py

Cases are grouped by laser power and scan speed. Within a group, an arrow
joins each mesh to the next finer one, and the finest mesh (or the only one) has a filled
marker. A case whose name ends in "_baseResistivity" is the case of the same
name before the calibration.
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
BASE = "_baseResistivity"
GREY = "#8a8a86"
PALETTE = ["#1f5fa8", "#c0462c", "#2c8a5a", "#7a4fa3", "#b8860b"]
ROLE = {(195.0, 1.2): "calibration: the resistivity was chosen at this speed",
        (195.0, 1.0): "validation: same resistivity, nothing adjusted"}


def rows(name):
    with open(HERE / "results" / name) as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


experiment = {(float(r["power_W"]), float(r["scan_speed_m_s"])):
              (float(r["width_um"]), float(r["depth_um"])) for r in rows("dilip2017.csv")}
cases = [dict(name=r["case"], key=(float(r["power"]), float(r["scan_speed"])),
              cell=1e6*float(r["cell_size"]), fill=100*float(r["fill_mean"]),
              width=1e6*float(r["width_mean"]), widthSd=1e6*float(r["width_sd"]),
              depth=1e6*float(r["depth_mean"]), depthSd=1e6*float(r["depth_sd"]))
         for r in rows("simtopcSummary.csv")]

print(f"{'case':<30}{'width (um)':>16}{'depth (um)':>16}{'fill':>7}{'experiment':>13}{'difference':>20}")
for c in cases:
    line = (f"{c['name']:<30}{c['width']:>9.1f} +- {c['widthSd']:.1f}{c['depth']:>10.1f} +- {c['depthSd']:.1f}"
            f"{c['fill']:>6.1f}%")
    if c["key"] in experiment:
        w, d = experiment[c["key"]]
        line += f"{w:>7.0f} / {d:.0f}{100*(c['width']/w - 1):>+11.1f} / {100*(c['depth']/d - 1):+.1f} %"
    print(line)

fig, ax = plt.subplots(figsize=(7.4, 5.8), dpi=160)
groups = sorted({c["key"] for c in cases}, reverse=True)
for n, key in enumerate(groups):
    col = PALETTE[n % len(PALETTE)]
    label = f"{key[1]:.1f} m/s" + (f", {ROLE[key]}" if key in ROLE else f", {key[0]:.0f} W")
    if key in experiment:
        w, d = experiment[key]
        ax.fill([0.9*w, 1.1*w, 1.1*w, 0.9*w], [0.9*d, 0.9*d, 1.1*d, 1.1*d], color=col, alpha=0.10, lw=0)
        ax.plot(w, d, "*", ms=17, color=col, mec="black", mew=0.8, zorder=4)
        ax.text(0.9*w + 0.6, 1.1*d - 0.8, label, color=col, fontsize=8.5, va="top")
    meshes = sorted((c for c in cases if c["key"] == key and not c["name"].endswith(BASE)),
                    key=lambda c: -c["cell"])
    for coarse, fine in zip(meshes, meshes[1:]):
        ax.annotate("", (fine["width"], fine["depth"]), (coarse["width"], coarse["depth"]),
                    arrowprops=dict(arrowstyle="->", color=col, lw=1.3, shrinkA=6, shrinkB=6))
    for c in meshes:
        finest = c is meshes[-1]
        ax.errorbar(c["width"], c["depth"], xerr=c["widthSd"], yerr=c["depthSd"], fmt="o", ms=8,
                    color=col, mfc=col if finest else "white", mec=col, mew=1.6, elinewidth=1,
                    capsize=2.5, zorder=3)
        if key not in experiment:
            ax.text(c["width"] + 1, c["depth"] + 1, label, color=col, fontsize=8.5)
for c in cases:
    if c["name"].endswith(BASE):
        after = [k for k in cases if k["name"] == c["name"][:-len(BASE)]]
        if after:
            ax.annotate("", (after[0]["width"], after[0]["depth"]), (c["width"], c["depth"]),
                        arrowprops=dict(arrowstyle="->", color=GREY, lw=1.3, ls="--", shrinkA=6, shrinkB=6))
        ax.errorbar(c["width"], c["depth"], xerr=c["widthSd"], yerr=c["depthSd"], fmt="o", ms=8,
                    color=GREY, elinewidth=1, capsize=2.5, zorder=3)
        ax.text(c["width"] - c["widthSd"] - 0.8, c["depth"], "before the calibration (1.9e-6 Ω m)",
                color="#5f5f5b", fontsize=8.5, ha="right", va="center")
h = [plt.Line2D([], [], marker="*", ms=13, ls="none", color="#666", mec="black", mew=0.8,
                label="experiment (Dilip et al. 2017)"),
     plt.Line2D([], [], marker="o", ms=8, ls="none", mfc="white", mec="#666", mew=1.6,
                label="laserbeamFoam, coarser mesh"),
     plt.Line2D([], [], marker="o", ms=8, ls="none", color="#666", label="laserbeamFoam, finest mesh")]
ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3, frameon=False, fontsize=8.5,
          handletextpad=0.3, columnspacing=1.4)
ax.set_ylim(bottom=min(ax.get_ylim()[0], min(c["depth"] for c in cases) - 3))
ax.set_xlabel("melt-pool width (µm)")
ax.set_ylabel("melt-pool depth (µm)")
ax.set_title("Ti-6Al-4V single tracks, 195 W: measured with SimToPC (mean ± sd along the track)\n"
             "Coloured points: calibrated electrical resistivity, 2.5e-6 Ω m; boxes: experiment ± 10 %",
             fontsize=9.5, loc="left")
ax.grid(color="#e8e8e6", lw=0.6)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(HERE / "validation.png")
print("written validation.png")
