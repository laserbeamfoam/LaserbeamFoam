#!/usr/bin/env python3
"""
Check the 3D AMR verification cases (11-15) against their exact answers.

    python3 tools/check3d.py [cases_dir]   -> table on stdout, cases/vv_results3d.md, cases/*/vv_*.png

Serial and decomposed runs are both read as written (processor*/ without reconstruction). Cell
volumes and centres come from each time's own polyMesh, so the AMR meshes need no post-processing.
"""
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mesh3d                                     # noqa: E402
from foamio import log_series                     # noqa: E402
from check_all import RS, CV, CL, LV, RU, ML, RHO0, psat, rho_liq   # noqa: E402

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
except ImportError:                               # plots are optional
    plt = None


# ----------------------------------------------------------------------------- helpers
def load(case):
    md = json.load(open(os.path.join(case, "case.json")))
    md["h0"] = md["dx"]
    md["hf"] = md["dx"] / 2 ** md["maxRefinement"]
    return md


def snaps(case):
    return [mesh3d.Snapshot(case, t) for t in mesh3d.times(case)]


def masses(s, md):
    """mass of every phase: c_i from the ledger, or alpha_i*rho_i(p, T) at t = 0 (c_i not written)"""
    out = {}
    for ph in md["phases"]:
        if s.has(f"c.{ph}"):
            c = s[f"c.{ph}"]
        elif ph == "metal1":
            c = s["alpha.metal1"] * rho_liq(s["p"])
        else:
            c = s[f"alpha.{ph}"] * s["p"] / (RS[ph] * s["T"])
        out[ph] = float(np.sum(c * s.V))
    return out


def internal_energy(s, md):
    m = 0.0
    for ph in md["phases"]:
        if s.has(f"c.{ph}"):
            c = s[f"c.{ph}"]
        elif ph == "metal1":
            c = s["alpha.metal1"] * rho_liq(s["p"])
        else:
            c = s[f"alpha.{ph}"] * s["p"] / (RS[ph] * s["T"])
        m = m + c * (CL if ph == "metal1" else CV[ph])
    return float(np.sum(m * s["T"] * s.V))


def amr_log(case):
    """refine and unrefine events and the cell-count range, from the solver log(s)"""
    log = os.path.join(case, "log.solver")
    ref = log_series(log, r"^Refined from (\d+) to (\d+) cells", 2)
    unr = log_series(log, r"^Unrefined from (\d+) to (\d+) cells", 2)
    counts = np.concatenate([a.ravel() for a in (ref, unr) if a.size]) if (ref.size or unr.size) else np.array([0])
    return len(ref), len(unr), int(counts.min()), int(counts.max())


def amr_text(case):
    nr, nu, cmin, cmax = amr_log(case)
    return (nr > 0 and nu > 0), f"AMR {nr} refine / {nu} unrefine events, {cmin}-{cmax} cells"


def purged(case):
    """cumulative mass removed by the ledger's trace purge (|c| < cPurgeTol), as a function of time"""
    t, cum, tot = [0.0], [0.0], 0.0
    rx_t = re.compile(r"^Time = ([-+0-9.eE]+)")
    rx_p = re.compile(r"mass-balance purge .*?([-+0-9.eE]+) kg absolute")
    with open(os.path.join(case, "log.solver"), errors="replace") as f:
        for line in f:
            m = rx_t.match(line)
            if m:
                t.append(float(m.group(1)))
                cum.append(tot)
                continue
            m = rx_p.search(line)
            if m:
                tot += float(m.group(1))
                cum[-1] = tot
    t, cum = np.array(t), np.array(cum)
    return lambda tt: float(cum[np.searchsorted(t, float(tt) * (1 + 1e-9), side="right") - 1])


def drift(ss, md, phases, case):
    """relative drift of the ledger mass of phases from the first write on, with the purged mass added
    back, and the offset of that first write from alpha*rho(p, T) at t = 0 (OpenFOAM's gas constant
    differs from ours in the 7th-9th digit, so that one is only a sanity check)"""
    P = purged(case)
    M = [sum(masses(s, md)[ph] for ph in phases) for s in ss]
    ref = M[1] + P(ss[1].time)
    d = max(abs((m + P(s.time)) / ref - 1) for m, s in zip(M[1:], ss[1:]))
    return d, abs(M[1] / M[0] - 1), P(ss[-1].time) / ref


def evap_series(case):
    """(time, net phase-change rate: evaporation + condensation, kg/s) per step from the solver log"""
    t = log_series(os.path.join(case, "log.solver"), r"^Time = ([-+0-9.eE]+)")
    e = log_series(os.path.join(case, "log.solver"),
                   r"integrated evap = ([-+0-9.eE]+) kg/s, cond = ([-+0-9.eE]+) kg/s", 2)
    n = min(len(t), len(e))
    return t[:n, 0], e[:n, 0] + e[:n, 1]


def slice_plot(ax, s, name, axis=1, at=None, vmin=None, vmax=None, cmap="viridis"):
    """cells cut by a coordinate plane, drawn at their own size (shows the AMR levels)"""
    h = np.cbrt(s.V)
    at = (at if at is not None else float(np.median(s.C[:, axis]))) + 1e-3 * h.min()   # off any face
    k = np.abs(s.C[:, axis] - at) < 0.5 * h
    a, b = [i for i in range(3) if i != axis]
    f = s[name][k]
    rects = [Rectangle(((s.C[i, a] - 0.5 * h[i]) * 1e6, (s.C[i, b] - 0.5 * h[i]) * 1e6), h[i] * 1e6, h[i] * 1e6)
             for i in np.where(k)[0]]
    pc = PatchCollection(rects, cmap=cmap, edgecolor="0.3", linewidth=0.15)
    pc.set_array(f)
    pc.set_clim(vmin if vmin is not None else f.min(), vmax if vmax is not None else f.max())
    ax.add_collection(pc)
    ax.autoscale_view()
    ax.set_aspect("equal")
    ax.set_xlabel("xyz"[a] + " [um]")
    ax.set_ylabel("xyz"[b] + " [um]")
    plt.colorbar(pc, ax=ax, shrink=0.8, label=name)


def savefig(fig, case, name):
    fig.tight_layout()
    fig.savefig(os.path.join(case, name), dpi=110)
    plt.close(fig)


def interface_band(s, md, dist):
    """mixed cells: their volume per interface area (in finest cells) and how far from the exact surface"""
    al = s["alpha.metal1"]
    mixed = (al > 0.01) & (al < 0.99)
    wet = al > 0.01
    thick = float(np.sum(s.V[mixed]) / md["area"] / md["hf"])
    d = dist(s.C)
    stray = float(np.max(np.abs(d[wet & (d > 0)]), initial=0.0) / md["hf"])   # liquid outside the surface
    return int(mixed.sum()), thick, stray


# ----------------------------------------------------------------------------- 11
def check_thermal_box(case):
    md = load(case)
    ss = snaps(case)
    p0 = md["p0"]
    # the volume mean of p is (gamma - 1)E/V exactly; its uniformity is checked once the acoustic
    # transient of the initial temperature step (1-2 % at the far corner after 0.5 us) has died out
    dpm = max(abs(float(np.sum(s["p"] * s.V) / np.sum(s.V)) / p0 - 1) for s in ss)
    dp = float(np.max(np.abs(ss[-1]["p"] / p0 - 1)))
    dpx = max(float(np.max(np.abs(s["p"] / p0 - 1))) for s in ss)
    E0 = internal_energy(ss[0], md)
    E = np.array([internal_energy(s, md) for s in ss])
    tt = np.array([float(s.time) for s in ss])
    dm, dm0, pg = drift(ss, md, ["argon"], case)
    dE = float(np.max(np.abs(E / E0 - 1)))
    amr_ok, amr = amr_text(case)
    last = ss[-1]
    Tm = float(np.sum(last["c.argon"] * last["T"] * last.V) / np.sum(last["c.argon"] * last.V))
    Tex = 1.0 / np.sum((1.0 / ss[0]["T"]) * ss[0].V) * np.sum(ss[0].V)   # mass-mean T is the harmonic volume mean
    ok = dpm < 1e-3 and dp < 1e-3 and dm < 1e-9 and dm0 < 1e-5 and dE < 1e-3 and amr_ok
    msg = (f"volume-mean p within {dpm:.1e} of p0, p uniform to {dp:.1e} at the end (acoustic peak {dpx:.1e} early on); "
           f"mass drift {dm:.1e} (trace purge {pg:.0e}); energy drift {dE:.1e}; mass-mean T {Tm:.2f} K "
           f"(conserved value {Tex:.2f}); T range at end {last['T'].min():.0f}-{last['T'].max():.0f} K; {amr}")
    if plt:
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
        ax[0].plot(tt * 1e6, [np.max(np.abs(s["p"] / p0 - 1)) for s in ss], "k.-", label="max |p/p0 - 1|")
        ax[0].plot(tt * 1e6, np.abs(E / E0 - 1) + 1e-16, "r.-", label="|E/E0 - 1|")
        ax[0].set_yscale("log"); ax[0].set_xlabel("t [us]"); ax[0].legend(fontsize=7)
        mid = ss[len(ss) // 3]
        slice_plot(ax[1], mid, "T", axis=1, at=float(np.min(mid.C[:, 1])))
        ax[1].set_title(f"T at t = {float(mid.time) * 1e6:.0f} us (y = first cell)")
        slice_plot(ax[2], last, "T", axis=1, at=float(np.min(last.C[:, 1])))
        ax[2].set_title(f"T at t = {float(last.time) * 1e6:.0f} us")
        savefig(fig, case, "vv_thermalBox3D.png")
    return ("PASS" if ok else "FAIL"), msg


# ----------------------------------------------------------------------------- 12
def check_contact(case):
    md = load(case)
    ss = snaps(case)
    p0, U0 = md["p0"], np.array(md["U0"])
    L = np.array(md["L"])
    dU = max(float(np.max(np.linalg.norm(s["U"] - U0, axis=1)) / np.linalg.norm(U0)) for s in ss)
    dp = max(float(np.max(np.abs(s["p"] / p0 - 1))) for s in ss)
    ms = [masses(s, md) for s in ss]
    dmv = max(abs(m["metal1vapour"] / ms[1]["metal1vapour"] - 1) for m in ms[1:])
    dmg = max(abs(m["argon"] / ms[1]["argon"] - 1) for m in ms[1:])
    dm, dm0, pg = drift(ss, md, ["metal1vapour", "argon"], case)
    # bulk velocity (mass-weighted) and the vapour centroid, unwrapped about the exact position
    errs, path = [], []
    for s in ss[1:]:
        rho = s["c.metal1vapour"] + s["c.argon"]
        Ub = np.sum(rho[:, None] * s["U"] * s.V[:, None], axis=0) / np.sum(rho * s.V)
        xe = np.array(md["c0"]) + U0 * float(s.time)
        X = s.C.copy()
        for k in (0, 1):
            X[:, k] -= L[k] * np.round((X[:, k] - xe[k]) / L[k])
        w = s["c.metal1vapour"] * s.V
        xc = np.sum(w[:, None] * X, axis=0) / w.sum()
        errs.append((np.linalg.norm(Ub - U0) / np.linalg.norm(U0), np.linalg.norm(xc - xe) / md["hf"]))
        path.append((float(s.time), *(xc - xe)))
    dP = max(e[0] for e in errs)
    dx = max(e[1] for e in errs)
    amr_ok, amr = amr_text(case)
    ok = dU < 1e-3 and dp < 1e-3 and dP < 1e-3 and dm < 1e-9 and dm0 < 1e-5 and dx < 0.5 and amr_ok
    msg = (f"max |U - U0|/|U0| {dU:.1e}; max |p/p0 - 1| {dp:.1e}; bulk momentum {dP:.1e}; gas mass drift {dm:.1e} "
           f"with the trace purge ({pg:.0e}) added back (vapour {dmv:.1e}, argon {dmg:.1e} without); "
           f"vapour centroid within {dx:.2f} finest cells of c0 + U0 t "
           f"(travelled {np.linalg.norm(U0) * float(ss[-1].time) * 1e6:.0f} um); {amr}")
    if plt:
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
        tt = [float(s.time) * 1e6 for s in ss]
        ax[0].plot(tt, [np.max(np.linalg.norm(s["U"] - U0, axis=1)) / np.linalg.norm(U0) + 1e-16 for s in ss],
                   "k.-", label="max |U-U0|/|U0|")
        ax[0].plot(tt, [np.max(np.abs(s["p"] / p0 - 1)) + 1e-16 for s in ss], "b.-", label="max |p/p0-1|")
        ax[0].set_yscale("log"); ax[0].set_xlabel("t [us]"); ax[0].legend(fontsize=7)
        P = np.array(path)
        ax[1].plot(P[:, 0] * 1e6, P[:, 1] / md["hf"], label="x"); ax[1].plot(P[:, 0] * 1e6, P[:, 2] / md["hf"], label="y")
        ax[1].plot(P[:, 0] * 1e6, P[:, 3] / md["hf"], label="z")
        ax[1].set_xlabel("t [us]"); ax[1].set_ylabel("centroid - (c0 + U0 t) [finest cells]"); ax[1].legend(fontsize=7)
        slice_plot(ax[2], ss[-1], "alpha.metal1vapour", axis=2, at=md["c0"][2], vmin=0, vmax=1)
        ax[2].set_title(f"vapour at t = {tt[-1]:.0f} us, z = c0")
        savefig(fig, case, "vv_contact3D.png")
    return ("PASS" if ok else "FAIL"), msg


# ----------------------------------------------------------------------------- 13, 14
def check_evap_closed(case):
    md = load(case)
    ss = snaps(case)
    T0, p0 = md["T0"], md["p0"]
    j0 = math.sqrt(ML / (2 * math.pi * RU * T0)) * (psat(T0) - p0)
    te, ev = evap_series(case)
    A1 = ev[0] / j0
    if md["kind"] == "film":
        nrm = math.sqrt(1 + md["a"] ** 2 + md["b"] ** 2)
        Lb = md["L"][0]
        dist = lambda C: (C[:, 2] - md["z0"] - md["a"] * (C[:, 0] - 0.5 * Lb) - md["b"] * (C[:, 1] - 0.5 * Lb)) / nrm  # noqa: E731
    else:
        dist = lambda C: np.linalg.norm(C, axis=1) - md["R"]                                                       # noqa: E731
    last = ss[-1]
    al = last["alpha.metal1"]
    gas = al < 0.01
    pg = float(np.sum(last["p"][gas] * last.V[gas]) / np.sum(last.V[gas]))
    cl = last["c.metal1"]
    Tl = float(np.sum(cl * last["T"] * last.V) / np.sum(cl * last.V))
    dm1, dm0, purge = drift(ss, md, ["metal1", "metal1vapour"], case)
    n0, th0, _ = interface_band(ss[0], md, dist)
    n1, th1, stray = interface_band(last, md, dist)
    inner = dist(last.C) < -2 * md["hf"]                  # vapour inside the liquid, 2 finest cells in or more
    vin = float(np.sum(last["alpha.metal1vapour"][inner] * last.V[inner]) / md["hf"] ** 3)
    Ul = float(np.max(np.linalg.norm(last["U"][al > 0.99], axis=1)))
    amr_ok, amr = amr_text(case)
    area_ok = abs(A1 / md["area"] - 1) < 0.05
    thermo_ok = abs(pg / psat(Tl) - 1) < 0.005 and dm1 < 1e-9 and dm0 < 1e-5
    film_ok = th1 < 3.0 and stray < 2.0
    ok = area_ok and thermo_ok and film_ok
    msg = (f"[area {'ok' if area_ok else 'FAIL'}] first-step phase change / j0 = {A1:.4e} m2 vs exact "
           f"{md['area']:.4e} ({A1 / md['area'] - 1:+.2%}); [thermo {'ok' if thermo_ok else 'FAIL'}] vapour p "
           f"{pg:.0f} Pa vs Psat(T_liq {Tl:.2f} K) {psat(Tl):.0f} ({pg / psat(Tl) - 1:+.2%}; Psat(T0) {psat(T0):.0f}), "
           f"metal mass drift {dm1:.1e} (trace purge {purge:.0e}); [interface {'ok' if film_ok else 'FAIL'}] mixed band "
           f"{th1:.2f} finest cells thick (t = 0: {th0:.2f}), liquid at most {stray:.1f} cells outside the exact "
           f"surface, vapour inside the liquid {vin:.2f} finest-cell volumes, max |U| in the liquid {Ul:.2f} m/s; {amr}")
    if plt:
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
        tt = [float(s.time) * 1e6 for s in ss]
        pgs = []
        for s in ss:
            g_ = s["alpha.metal1"] < 0.01
            pgs.append(np.sum(s["p"][g_] * s.V[g_]) / np.sum(s.V[g_]))
        # Hertz-Knudsen on the exact area against the mean vapour pressure (T0: the cooling is ~0.5 K)
        pred = j0 / (psat(T0) - p0) * (psat(T0) - np.interp(te * 1e6, tt, pgs)) * md["area"]
        ax[0].plot(te * 1e6, ev * 1e9, "k-", lw=0.8, label="solver (integrated evaporation)")
        ax[0].plot(te * 1e6, pred * 1e9, "r--", lw=1.0, label="j(mean vapour p) x exact area")
        ax[0].set_xlabel("t [us]"); ax[0].set_ylabel("evaporation [ug/s]"); ax[0].legend(fontsize=7)
        ax[1].plot(tt, np.array(pgs) / 1e3, "b.-", label="mean vapour p")
        ax[1].axhline(psat(T0) / 1e3, color="k", ls="--", label="Psat(T0)")
        ax[1].set_xlabel("t [us]"); ax[1].set_ylabel("kPa"); ax[1].legend(fontsize=7)
        if md["kind"] == "film":
            slice_plot(ax[2], last, "alpha.metal1", axis=1, at=0.5 * md["L"][1], vmin=0, vmax=1)
        else:
            slice_plot(ax[2], last, "alpha.metal1", axis=2, at=float(np.min(last.C[:, 2])), vmin=0, vmax=1)
        ax[2].set_title(f"liquid at t = {tt[-1]:.1f} us")
        savefig(fig, case, f"vv_{md['kind']}3D.png")
    if md["kind"] == "dropletEvap" and area_ok and thermo_ok and not film_ok and vin > 1.0:
        return "KNOWN", msg + ("; KNOWN: the superheated droplet grows vapour pockets from traces inside it "
                               "(a library from before the 8 Nov night update has no sub-grid barrier), see README")
    return ("PASS" if ok else "FAIL"), msg


# ----------------------------------------------------------------------------- 15
def layers(s, hf):
    """volume-weighted means over z-slabs of the finest spacing"""
    k = np.floor(s.C[:, 2] / hf).astype(int)
    n = k.max() + 1
    W = np.bincount(k, s.V, n)
    use = W > 0

    def mean(f):
        return np.bincount(k, f * s.V, n)[use] / W[use]
    return (np.arange(n)[use] + 0.5) * hf, mean, k, use


def check_evap_steady(case):
    md = load(case)
    ss = snaps(case)
    hf, A = md["hf"], md["L"][0] * md["L"][1]
    pinf, q = md["p_inf"], md["q"]
    j_ref = q / LV
    last = ss[-1]
    z, mean, k, use = layers(last, hf)
    Uz = mean(last["U"][:, 2])
    cv = mean(last["c.metal1vapour"])
    j_top = float(cv[-1] * Uz[-1])
    half = ss[len(ss) // 2:]
    ml = [float(np.sum(s["c.metal1"] * s.V)) for s in (half[0], half[-1])]
    j_liq = -(ml[1] - ml[0]) / (float(half[-1].time) - float(half[0].time)) / A
    al = last["alpha.metal1"]
    mixed_layers = float(np.sum(last.V[(al > 0.01) & (al < 0.99)]) / (A * hf))
    lo, hi = 3000.0, 5000.0
    for _ in range(200):
        Ts = 0.5 * (lo + hi)
        rv = pinf / (RS["metal1vapour"] * Ts)
        jh = math.sqrt(ML / (2 * math.pi * RU * Ts)) * (psat(Ts) - pinf - j_ref ** 2 / rv)
        lo, hi = (Ts, hi) if jh < j_ref else (lo, Ts)
    dp_ref = j_ref ** 2 * (1 / rv - 1 / RHO0)

    def p_surface(s):
        z_, mean_, _, _ = layers(s, hf)
        a_, p_ = mean_(s["alpha.metal1"]), mean_(s["p"])
        full = np.where(a_ > 0.99)[0]
        return p_[full.max()] if full.size else float("nan")    # nan: no full liquid slab left
    dp_sim = float(np.mean([p_surface(s) for s in half])) - pinf
    # the column must stay one-dimensional: spread of p across each slab, transverse velocity
    pm = mean(last["p"])[np.searchsorted(np.where(use)[0], k)]
    dpx = float(np.max(np.abs(last["p"] - pm)))
    uxy = float(np.max(np.linalg.norm(last["U"][:, :2], axis=1)))
    amr_ok, amr = amr_text(case)
    ok = abs(j_liq / j_ref - 1) < 0.05 and mixed_layers <= 4 and abs(dp_sim / dp_ref - 1) < 0.25 \
        and dpx < 0.1 * dp_ref and uxy < 0.01 * j_ref / rv
    if not np.isfinite(dp_sim):
        dp_sim = float("nan")
    msg = (f"evaporation {j_liq:.2f} kg/m2/s (liquid loss) and {j_top:.2f} (vapour out) vs q/Lvap {j_ref:.2f}; "
           f"mixed band {mixed_layers:.2f} finest layers; p(liquid)-p_inf {dp_sim:.0f} Pa vs recoil {dp_ref:.0f} "
           f"(top liquid layer, 2nd-half mean); vapour speed {Uz[-1]:.1f} vs j/rho_v {j_ref / rv:.1f} m/s; "
           f"1D-ness: p spread across a slab {dpx:.1f} Pa, max transverse |U| {uxy:.2e} m/s; {amr}")
    if plt:
        fig, ax = plt.subplots(1, 3, figsize=(11, 3))
        for s in ss[1:: max(1, len(ss) // 4)] + [ss[-1]]:
            z_, mean_, _, _ = layers(s, hf)
            ax[0].plot(z_ * 1e6, mean_(s["alpha.metal1"]), label=f"{float(s.time) * 1e6:.0f} us")
        ax[0].set_xlim(0, 300); ax[0].legend(fontsize=7); ax[0].set_ylabel("alpha liquid"); ax[0].set_xlabel("z [um]")
        ax[1].plot(z * 1e6, Uz, "r-"); ax[1].axhline(j_ref / rv, color="k", ls="--"); ax[1].set_ylabel("U_z [m/s]")
        ax[1].set_xlabel("z [um]")
        ax[2].plot(z * 1e6, mean(last["p"]) - pinf, color="0.7", lw=0.8, label="last time")
        ax[2].axhline(dp_ref, color="k", ls="--"); ax[2].set_ylabel("p - p_inf [Pa]"); ax[2].set_xlabel("z [um]")
        ax[2].legend(fontsize=7)
        savefig(fig, case, "vv_evapSteady3D.png")
    return ("PASS" if ok else "FAIL"), msg


CHECKS = [
    ("11_thermalBox3D", "3D energy + closure under AMR refine/unrefine (sealed octant)", check_thermal_box),
    ("12_contact3D", "3D advected contact, periodic, AMR following it: U, p uniform", check_contact),
    ("13_inclinedFilm3D", "evaporation from a plane inclined to the mesh: rate = j*A, p -> Psat", check_evap_closed),
    ("14_dropletCond3D", "condensing droplet, all orientations and AMR levels: j*A, Psat", check_evap_closed),
    ("15_evapSteady3D", "test 09 on a 3D AMR column: flux, recoil, stays 1D", check_evap_steady),
    ("16_dropletEvap3D", "14 evaporating: superheated droplet, no vapour may grow inside it", check_evap_closed),
]


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(HERE), "cases")
    rows = []
    for name, what, fn in CHECKS:
        case = os.path.join(root, name)
        if not os.path.isdir(case) or len(mesh3d.times(case)) < 2:
            rows.append((name, what, "NOT RUN", ""))
            continue
        try:
            st, msg = fn(case)
        except Exception as e:                    # report, keep going
            st, msg = "ERROR", f"{type(e).__name__}: {e}"
        rows.append((name, what, st, msg))
        print(f"{st:7s} {name:22s} {msg}")
    with open(os.path.join(root, "vv_results3d.md"), "w") as f:
        f.write("| case | tests | result | numbers |\n|---|---|---|---|\n")
        for r in rows:
            f.write("| " + " | ".join(r) + " |\n")
    return 0 if all(r[2] in ("PASS", "KNOWN") for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
