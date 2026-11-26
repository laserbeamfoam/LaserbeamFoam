#!/usr/bin/env python3
"""
Check every verification case against its exact / semi-analytical answer.

    python3 tools/check_all.py [cases_dir]        -> table on stdout, vv_results.md, cases/*/vv_*.png

Each check returns PASS / FAIL / KNOWN (a deviation that is an identified model limitation,
documented in the README, with the number it should show) plus the measured numbers.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from foamio import times, read_field          # noqa: E402
from riemann import star, sample, shock_speed   # noqa: E402

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:                               # plots are optional
    plt = None

RR = 8314.47          # J/kmol/K, OpenFOAM's universal gas constant (specie R = RR/W)
RU = 8.314            # the solver's own constant in Clausius-Clapeyron / Hertz-Knudsen

# material constants exactly as gen_cases.py writes them
W = {"metal1vapour": 55.85, "argon": 39.948, "air": 28.96}
CP = {"metal1vapour": 2.5 * 8.314 / 0.05585, "argon": 2.5 * 8.314 / 0.039948, "air": 1005.0}
RS = {k: RR / w for k, w in W.items()}
CV = {k: CP[k] - RS[k] for k in W}
CL = 800.0            # liquid Fe Cv = Cp (adiabaticPerfectFluid)
LV, TB, ML = 6.09e6, 3134.0, 0.05585
RHO0, BL, GL = 7030.0, 1e7, 1.2


def psat(T):
    return 1e5 * np.exp(ML * LV / (TB * RU) * (1 - TB / T))


def rho_liq(p, B=BL):
    return RHO0 * ((p + B) / (1e5 + B)) ** (1 / GL)


def load(case):
    md = json.load(open(os.path.join(case, "case.json")))
    md["x"] = (np.arange(md["N"]) + 0.5) * md["dx"]
    return md


def F(case, t, name, md):
    v = read_field(case, t, name, md["N"])
    return v[:, 0] if v.ndim == 2 else v


def savefig(fig, case, name):
    if fig is not None:
        fig.tight_layout()
        fig.savefig(os.path.join(case, name), dpi=110)
        plt.close(fig)


# ----------------------------------------------------------------------------- checks
def check_acoustic(case):
    md = load(case)
    g = CP["argon"] / CV["argon"]
    cad, ciso = math.sqrt(g * RS["argon"] * md["T0"]), math.sqrt(RS["argon"] * md["T0"])
    ts, xs, ratio = [], [], []
    for t in times(case)[2:]:                     # skip t=0 and the still-overlapping first write
        p, T = F(case, t, "p", md), F(case, t, "T", md)
        right = md["x"] > md["x0"]
        i = np.where(right)[0][np.argmax(p[right])]
        y0, y1, y2 = p[i - 1], p[i], p[i + 1]
        xs.append(md["x"][i] + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2) * md["dx"])
        ts.append(float(t))
        ratio.append(((T[i] - md["T0"]) / md["T0"]) / ((p[i] - md["p0"]) / md["p0"]))
    c = np.polyfit(ts, xs, 1)[0]
    rth = np.mean(ratio)
    ok = abs(c / cad - 1) < 0.01 and abs(rth / ((g - 1) / g) - 1) < 0.01
    msg = (f"wave speed {c:.1f} m/s = {c / cad:.4f} c_adiabatic ({cad:.1f}; isothermal would be {ciso:.1f}); "
           f"(dT/T)/(dp/p) = {rth:.4f} vs (g-1)/g = {(g - 1) / g:.4f}")
    if plt:
        fig, ax = plt.subplots(figsize=(7, 3.2))
        for t in times(case)[::3]:
            ax.plot(md["x"] * 1e3, (F(case, t, "p", md) / md["p0"] - 1) * 100, label=f"t={float(t) * 1e6:.2f} us")
        ax.set_xlabel("x [mm]"); ax.set_ylabel("p/p0 - 1 [%]"); ax.legend(fontsize=7)
        savefig(fig, case, "vv_acoustic.png")
    return ("PASS" if ok else "FAIL"), msg


def check_thermal_box(case):
    md = load(case)
    tl = times(case)
    t0, tend = tl[1], tl[-1]
    c0, c1 = F(case, t0, "c.argon", md), F(case, tend, "c.argon", md)
    m0, m1 = c0.sum() * md["dx"], c1.sum() * md["dx"]
    T, p = F(case, tend, "p", md), None
    p = F(case, tend, "p", md)
    T = F(case, tend, "T", md)
    Tex = 2 / (1 / md["Thot"] + 1 / md["Tcold"])
    E0 = CV["argon"] * (md["p0"] / RS["argon"]) * md["L"]      # sum rho Cv T dx at t=0 (rho T = p/R in every cell)
    E1 = (c1 * CV["argon"] * T).sum() * md["dx"]
    pmax_dev = max(np.abs(F(case, t, "p", md) / md["p0"] - 1).max() for t in tl[1:])
    Tm = (c1 * T).sum() / c1.sum()
    ok = abs(m1 / m0 - 1) < 1e-9 and pmax_dev < 1e-3 and abs(E1 / E0 - 1) < 1e-3 and abs(Tm / Tex - 1) < 5e-3
    msg = (f"max|p/p0-1| over the run {pmax_dev:.1e} (exact 0); end T {T.min():.1f}-{T.max():.1f}, "
           f"mass-mean {Tm:.2f} K vs {Tex:.2f}; energy drift {E1 / E0 - 1:.1e}; mass drift {m1 / m0 - 1:.1e}")
    return ("PASS" if ok else "FAIL"), msg


def check_two_gas_box(case):
    md = load(case)
    tl = times(case)
    hot, cold = md["hot"], md["cold"]
    ch, cc = F(case, tl[1], "c." + hot, md), F(case, tl[1], "c." + cold, md)
    mh, mc = ch.sum() * md["dx"], cc.sum() * md["dx"]
    E0 = mh * CV[hot] * md["Thot"] + mc * CV[cold] * md["Tcold"]
    Tf = E0 / (mh * CV[hot] + mc * CV[cold])
    pf = (mh * RS[hot] + mc * RS[cold]) * Tf / md["L"]
    p, T = F(case, tl[-1], "p", md), F(case, tl[-1], "T", md)
    ch1, cc1 = F(case, tl[-1], "c." + hot, md), F(case, tl[-1], "c." + cold, md)
    E1 = ((ch1 * CV[hot] + cc1 * CV[cold]) * T).sum() * md["dx"]
    ok = abs(p.mean() / pf - 1) < 1e-3 and abs(T.mean() / Tf - 1) < 5e-3 and abs(E1 / E0 - 1) < 1e-3
    msg = (f"end p {p.mean():.1f} Pa vs exact {pf:.1f} ({p.mean() / pf - 1:+.1e}); end T {T.min():.2f}-{T.max():.2f} "
           f"vs {Tf:.2f}; energy drift {E1 / E0 - 1:.1e} (p0 {md['p0']:.0f}: mixing a hot monatomic gas into a "
           f"cold diatomic one at fixed volume lowers p)")
    if plt:
        fig, ax = plt.subplots(1, 2, figsize=(8, 3))
        tt = [float(t) for t in tl[1:]]
        ax[0].plot(np.array(tt) * 1e3, [F(case, t, "p", md).mean() for t in tl[1:]], "o-")
        ax[0].axhline(pf, color="k", ls="--"); ax[0].set_xlabel("t [ms]"); ax[0].set_ylabel("p [Pa]")
        for t in tl[::3]:
            ax[1].plot(md["x"] * 1e3, F(case, t, "T", md), label=f"{float(t) * 1e3:.1f} ms")
        ax[1].axhline(Tf, color="k", ls="--"); ax[1].set_xlabel("x [mm]"); ax[1].set_ylabel("T [K]")
        ax[1].legend(fontsize=7)
        savefig(fig, case, "vv_twoGasBox.png")
    return ("PASS" if ok else "FAIL"), msg


def check_sod(case):
    md = load(case)
    g = CP["argon"] / CV["argon"]
    R = RS["argon"]
    rL, rR = md["pL"] / (R * md["TL"]), md["pR"] / (R * md["TR"])
    ps, us = star(rL, 0, md["pL"], rR, 0, md["pR"], g)
    t = times(case)[-1]
    tt = float(t)
    x = md["x"]
    p, U, T = F(case, t, "p", md), F(case, t, "U", md), F(case, t, "T", md)
    sR = shock_speed(rR, 0, md["pR"], ps, g)
    xc, xs = md["x0"] + us * tt, md["x0"] + sR * tt
    mid = (x > xc + 0.25 * (xs - xc)) & (x < xs - 0.25 * (xs - xc))
    r2 = rR * (ps / md["pR"] + (g - 1) / (g + 1)) / ((g - 1) / (g + 1) * ps / md["pR"] + 1)
    r2i = rR * (ps / md["pR"]) ** (1 / g)
    half = 0.5 * (ps + md["pR"])
    j = np.where((x > xc) & (p < half))[0][0]
    xsim = x[j - 1] + (half - p[j - 1]) / (p[j] - p[j - 1]) * md["dx"]
    eps, epu = p[mid].mean() / ps - 1, U[mid].mean() / us - 1
    esh = (xsim - md["x0"]) / (xs - md["x0"]) - 1
    ok = abs(eps) < 0.02 and abs(epu) < 0.02 and abs(esh) < 0.02
    TRH, Tisen = ps / (R * r2), ps / (R * r2i)
    msg = (f"p* {p[mid].mean():.0f} vs {ps:.0f} ({eps:+.1%}), u* {U[mid].mean():.1f} vs {us:.1f} ({epu:+.1%}), "
           f"shock position {esh:+.1%}; post-shock T {T[mid].mean():.1f} K vs Rankine-Hugoniot {TRH:.1f} "
           f"(an isentropic compression would give {Tisen:.1f})")
    if plt:
        ref = sample(x, tt, md["x0"], rL, 0, md["pL"], rR, 0, md["pR"], g)
        fig, ax = plt.subplots(1, 3, figsize=(11, 3))
        for a, (sim, r, lab) in zip(ax, [(p, ref[2], "p [Pa]"), (U, ref[1], "u [m/s]"), (T, ref[2] / (R * ref[0]), "T [K]")]):
            a.plot(x * 1e3, r, "k-", lw=1, label="exact"); a.plot(x * 1e3, sim, "r.", ms=2, label="solver")
            a.set_xlabel("x [mm]"); a.set_ylabel(lab)
        ax[0].legend(fontsize=7)
        savefig(fig, case, "vv_sod.png")
    return ("PASS" if ok else "FAIL"), msg


def check_contact(case, air):
    md = load(case)
    tl = times(case)
    devU = max(np.abs(F(case, t, "U", md) - md["U0"]).max() / md["U0"] for t in tl[1:])
    devp = max(np.abs(F(case, t, "p", md) / md["p0"] - 1).max() for t in tl[1:])
    gas = md["gas"]
    m = lambda t: sum(F(case, t, "c." + ph, md) for ph in ("metal1vapour", gas))
    mom0 = (m(tl[1]) * F(case, tl[1], "U", md)).sum()
    mom1 = (m(tl[-1]) * F(case, tl[-1], "U", md)).sum()
    if not air:
        ok = devp < 1e-3 and devU < 0.05
        status = "PASS" if (devp < 1e-3 and devU < 1e-3) else ("KNOWN" if ok else "FAIL")
        msg = (f"max|p/p0-1| {devp:.1e}, max|U/U0-1| {devU:.1e} (exact 0 for equal gamma); "
               f"momentum change over the run {mom1 / mom0 - 1:+.2%}")
        return status, msg
    # full-mixing bound for the pressure (constant volume, energy conserved)
    T0 = F(case, "0", "T", md)
    av = F(case, "0", "alpha.metal1vapour", md)
    nv = (av * md["p0"] / (RR * T0)).sum()
    na = ((1 - av) * md["p0"] / (RR * T0)).sum()
    E = (av * md["p0"] / (CP["metal1vapour"] / CV["metal1vapour"] - 1)).sum() + \
        ((1 - av) * md["p0"] / (CP[gas] / CV[gas] - 1)).sum()
    C = nv * CV["metal1vapour"] * W["metal1vapour"] + na * CV[gas] * W[gas]
    pmix = (nv + na) * RR * (E / C) / md["N"]
    pend = F(case, tl[-1], "p", md).mean()
    ok = pmix * (1 - 1e-3) <= pend <= md["p0"] * (1 + 1e-3) and devU < 0.1
    msg = (f"max|U/U0-1| {devU:.1e}; mean p fell to {pend:.0f} Pa (full-mixing limit {pmix:.0f}, p0 {md['p0']:.0f}): "
           f"cells mixed by upwind diffusion contract, as the exact single-T mixing of different gammas does; "
           f"momentum change {mom1 / mom0 - 1:+.2%}")
    return ("KNOWN" if ok else "FAIL"), msg


def film_state(case, md, t, h0):
    """number of mixed cells and how far liquid has been carried above the initial film"""
    al = F(case, t, "alpha.metal1", md)
    mixed = int(((al > 0.01) & (al < 0.99)).sum())
    wet = md["x"][al > 0.01]
    spread = (wet.max() - h0) / md["dx"] if len(wet) else 0.0
    return mixed, spread


def gas_pressure(case, md, t):
    al = F(case, t, "alpha.metal1", md)
    p = F(case, t, "p", md)
    gas = al < 0.01
    return float(np.median(p[gas])) if gas.any() else float("nan"), float(p.max())


def check_evap_iso(case):
    md = load(case)
    tl = times(case)
    pg, pmax = gas_pressure(case, md, tl[-1])
    cl0, cv0 = F(case, tl[1], "c.metal1", md), F(case, tl[1], "c.metal1vapour", md)
    cl, cv = F(case, tl[-1], "c.metal1", md), F(case, tl[-1], "c.metal1vapour", md)
    mtot0, mtot = (cl0 + cv0).sum(), (cl + cv).sum()
    ps = psat(md["T0"])
    mixed, spread = film_state(case, md, tl[-1], md["h_l"])
    thermo = abs(pg / ps - 1) < 0.01 and abs(mtot / mtot0 - 1) < 1e-9
    film = mixed <= 4 and spread < 3
    msg = (f"[thermo {'ok' if thermo else 'FAIL'}] vapour p {pg:.0f} Pa vs Psat(T0) {ps:.0f} ({pg / ps - 1:+.2%}), "
           f"metal mass drift {mtot / mtot0 - 1:.1e}; [film {'ok' if film else 'FAIL'}] {mixed} mixed cells, liquid "
           f"carried {spread:.0f} cells above the film, max p {pmax:.0f} Pa (in trapped liquid)")
    if plt:
        plot_film(case, md, "vv_evapBox.png")
    return ("PASS" if thermo and film else "FAIL"), msg


def plot_film(case, md, name):
    fig, ax = plt.subplots(1, 2, figsize=(9, 3))
    tl = times(case)
    for t in tl[:: max(1, len(tl) // 5)] + [tl[-1]]:
        ax[0].plot(md["x"] * 1e6, F(case, t, "alpha.metal1", md), label=f"{float(t) * 1e6:.0f} us")
        ax[1].plot(md["x"] * 1e6, F(case, t, "p", md) / 1e3)
    ax[0].set_xlabel("x [um]"); ax[0].set_ylabel("alpha liquid"); ax[0].legend(fontsize=7)
    ax[1].set_xlabel("x [um]"); ax[1].set_ylabel("p [kPa]")
    savefig(fig, case, name)


def evap_adiabatic_reference(md, ml0, mv0):
    """integrate dT/dm_v = -LvapE/C_tot (sealed, adiabatic, uniform p) until p(T, m_v) = Psat(T)"""
    L = md["L"]
    T, mv, ml = md["T0"], mv0, ml0
    dm = 1e-9
    for _ in range(10 ** 7):
        Vg = L - ml / rho_liq(1e5)
        p = mv * RS["metal1vapour"] * T / Vg
        if p >= psat(T):
            return T, p, mv
        LvE = LV - p * (Vg / mv - 1 / rho_liq(p))
        C = ml * CL + mv * CV["metal1vapour"]
        T -= LvE * dm / C
        mv += dm
        ml -= dm
    raise RuntimeError("no equilibrium")


def check_evap_adiabatic(case):
    md = load(case)
    tl = times(case)
    ml0 = (F(case, "0", "alpha.metal1", md) * rho_liq(md["p0"])).sum() * md["dx"]
    av0 = F(case, "0", "alpha.metal1vapour", md)
    mv0 = (av0 * md["p0"] / (RS["metal1vapour"] * md["T0"])).sum() * md["dx"]
    Tf, pf, mvf = evap_adiabatic_reference(md, ml0, mv0)
    pg, pmax = gas_pressure(case, md, tl[-1])
    T = F(case, tl[-1], "T", md)
    cv = F(case, tl[-1], "c.metal1vapour", md)
    mixed, spread = film_state(case, md, tl[-1], md["h_l"])
    dT_ref = md["T0"] - Tf
    thermo = abs(pg / pf - 1) < 0.005 and abs((md["T0"] - np.median(T)) / dT_ref - 1) < 0.25
    film = mixed <= 4 and spread < 3
    msg = (f"[thermo {'ok' if thermo else 'FAIL'}] vapour p {pg:.0f} Pa vs Psat(T_f) {pf:.0f} ({pg / pf - 1:+.2%}; "
           f"Psat(T0) would be {psat(md['T0']):.0f}); cooling {md['T0'] - np.median(T):.2f} K vs {dT_ref:.2f} K; "
           f"vapour mass {cv.sum() * md['dx']:.4e} vs {mvf:.4e}; [film {'ok' if film else 'FAIL'}] {mixed} mixed cells, "
           f"liquid carried {spread:.0f} cells up")
    if plt:
        plot_film(case, md, "vv_evapBox.png")
    return ("PASS" if thermo and film else "FAIL"), msg


def check_evap_steady(case):
    md = load(case)
    tl = times(case)
    t = tl[-1]
    x = md["x"]
    al, U, p, T = F(case, t, "alpha.metal1", md), F(case, t, "U", md), F(case, t, "p", md), F(case, t, "T", md)
    cv = F(case, t, "c.metal1vapour", md)
    q = md["q"]
    j_ref = q / LV
    j_top = cv[-1] * U[-1]
    # liquid mass loss rate over the last half of the run
    t1, t2 = tl[len(tl) // 2], tl[-1]
    ml = lambda tt: F(case, tt, "c.metal1", md).sum() * md["dx"]
    j_liq = -(ml(t2) - ml(t1)) / (float(t2) - float(t1))
    mixed = int(((al > 0.01) & (al < 0.99)).sum())
    pinf = md["p_inf"]
    # reference surface temperature: Hertz-Knudsen gives the reference flux against the vapour-side
    # pressure p_inf + j^2/rho_v (pure vapour leaving the surface from rest)
    lo, hi = 3000.0, 5000.0
    for _ in range(200):
        Ts = 0.5 * (lo + hi)
        rv = pinf / (RS["metal1vapour"] * Ts)
        jh = math.sqrt(ML / (2 * math.pi * RU * Ts)) * (psat(Ts) - pinf - j_ref ** 2 / rv)
        lo, hi = (Ts, hi) if jh < j_ref else (lo, Ts)
    dp_ref = j_ref ** 2 * (1 / rv - 1 / RHO0)
    # liquid pressure at the interface (the jump condition holds there): the topmost full liquid
    # cell, averaged over the second half of the run. The soft liquid (B = 1e7 Pa, c = 41 m/s)
    # rings at its quarter-wave period 4h/c = 10 us with the pressure node at the free surface
    def p_surface(tt):
        a_, p_ = F(case, tt, "alpha.metal1", md), F(case, tt, "p", md)
        return p_[np.where(a_ > 0.99)[0].max()]
    dp_sim = float(np.mean([p_surface(tt) for tt in tl[len(tl) // 2:]])) - pinf
    ok = abs(j_liq / j_ref - 1) < 0.05 and mixed <= 4 and abs(dp_sim / dp_ref - 1) < 0.25
    msg = (f"evaporation {j_liq:.2f} kg/m2/s (liquid loss) and {j_top:.2f} (vapour out) vs q/Lvap {j_ref:.2f}; "
           f"{mixed} mixed cells (a sharp film has 1-3); p(liquid)-p_inf {dp_sim:.0f} Pa vs recoil j^2/rho_v "
           f"{dp_ref:.0f} (top liquid cell, 2nd-half mean); vapour speed {U[-1]:.1f} vs j/rho_v {j_ref / rv:.1f} m/s; "
           f"T_s ref (Hertz-Knudsen) {Ts:.1f} K")
    if plt:
        fig, ax = plt.subplots(1, 3, figsize=(11, 3))
        for tt in tl[1:: max(1, len(tl) // 4)] + [tl[-1]]:
            ax[0].plot(x * 1e6, F(case, tt, "alpha.metal1", md), label=f"{float(tt) * 1e6:.0f} us")
        ax[0].legend(fontsize=7); ax[0].set_ylabel("alpha liquid"); ax[0].set_xlabel("x [um]")
        ax[1].plot(x * 1e6, U, "r-"); ax[1].axhline(j_ref / rv, color="k", ls="--"); ax[1].set_ylabel("U [m/s]")
        ax[1].set_xlabel("x [um]")
        pm = np.mean([F(case, tt, "p", md) for tt in tl[len(tl) // 2:]], axis=0)
        ax[2].plot(x * 1e6, p - pinf, color="0.7", lw=0.8, label="last time")
        ax[2].plot(x * 1e6, pm - pinf, "g-", label="2nd-half mean")
        ax[2].axhline(dp_ref, color="k", ls="--"); ax[2].set_ylabel("p - p_inf [Pa]"); ax[2].legend(fontsize=7)
        ax[2].set_xlabel("x [um]")
        savefig(fig, case, "vv_evapSteady.png")
    return ("PASS" if ok else "FAIL"), msg


def check_stefan(case):
    from stefan_reference import Stefan
    md = load(case)
    S = Stefan(rho=md["rho"], k=md["k"], c=md["cp"], L=md["L"], Ts=md["Ts"], Tl=md["Tl"], Tw=md["Tw"], Ti=md["Ti"])
    worst_f, worst_T = 0.0, 0.0
    x = md["x"]
    for t in times(case)[1:]:
        tt = float(t)
        e, T = F(case, t, "epsilon1", md), F(case, t, "T", md)
        i = np.where(e < 0.5)[0][0]
        xs = x[i - 1] + (0.5 - e[i - 1]) / (e[i] - e[i - 1]) * md["dx"]
        xr = S.fronts(tt)[2]
        worst_f = max(worst_f, abs(xs / xr - 1))
        far = x < 0.6 * md["L"]
        worst_T = max(worst_T, np.sqrt(np.mean((T[far] - S.T(x[far], tt)) ** 2)) / abs(md["Tw"] - md["Ti"]))
    ok = worst_f < 0.01 and worst_T < 0.01
    msg = f"eps=0.5 front within {worst_f:.2%} of the similarity solution, T L2 error {worst_T:.3%} of (Tw-Ti)"
    return ("PASS" if ok else "FAIL"), msg


# ----------------------------------------------------------------------------- liquid-vapour Stefan problems
def liquid_mass(case, t, md):
    """liquid metal per unit area [kg/m2] from the ledger (at t = 0, before the solver writes one,
    from alpha and the liquid equation of state, which is how the ledger starts)"""
    if float(t) == 0:
        return float((F(case, t, "alpha.metal1", md) * rho_liq(F(case, t, "p", md))).sum() * md["dx"])
    return float(F(case, t, "c.metal1", md).sum() * md["dx"])


def vapour_mass(case, t, md):
    """metal vapour per unit area [kg/m2] from the ledger (at t = 0 from alpha and the perfect gas)"""
    if float(t) == 0:
        rv = F(case, t, "p", md) * W["metal1vapour"] / (RR * F(case, t, "T", md))
        return float((F(case, t, "alpha.metal1vapour", md) * rv).sum() * md["dx"])
    return float(F(case, t, "c.metal1vapour", md).sum() * md["dx"])


def first_crossing(x, a, level=0.5):
    """position where a first rises through level, scanning from x[0]"""
    i = int(np.argmax(a >= level))
    if i == 0:
        return float("nan")
    return float(x[i - 1] + (level - a[i - 1]) / (a[i] - a[i - 1]) * (x[i] - x[i - 1]))


def liquid_cells(case, t, md):
    """liquid metal per cell [kg/m3] from the ledger (at t = 0 from alpha and the equation of state)"""
    if float(t) == 0:
        return F(case, t, "alpha.metal1", md) * rho_liq(F(case, t, "p", md))
    return F(case, t, "c.metal1", md)


def check_stefan_boiling(case):
    """Film boiling against Welch & Wilson's similarity solution: the film's growth and temperature,
    and where the heat the film conducts to the slab goes.

    Hertz-Knudsen needs the slab's surface 0.1-0.4 K above Tsat to evaporate at this rate, and the
    interface cell's single temperature hands that superheat to the liquid in it as the slab moves
    through the mesh: the slab warms by ~0.05 K and takes ~3 % of the heat (the same on 400 cells, at
    dt/4 and with a 100x less conductive slab). The reference's liquid stays at Tsat, so the evaporated
    mass falls short by that 3 % and the film lags by half of it; the energy check counts both."""
    from stefan_lv import FilmBoiling
    md = load(case)
    S = FilmBoiling(k=md["k_v"], rho=md["rho_v"], cp=md["cp_v"], Tw=md["Tw"], Tsat=md["Tsat"])
    tl = times(case)
    x, dx = md["x"], md["dx"]
    cl, Lv = md.get("c_l", 800.0), md.get("L_v", 6.09e6)

    def liquid_heat(t):
        """the liquid's sensible heat above Tsat [J/m2] (Cv = Cp for adiabaticPerfectFluid)"""
        return float(np.sum(liquid_cells(case, t, md) * cl * (F(case, t, "T", md) - md["Tsat"])) * dx)

    m0, h0 = liquid_mass(case, tl[0], md), liquid_heat(tl[0])
    worst_d, worst_T, worst_e, rows = 0.0, 0.0, 0.0, []
    for t in tl[1:]:
        tt = md["t0"] + float(t)
        al, T = F(case, t, "alpha.metal1", md), F(case, t, "T", md)
        xi = first_crossing(x, al)
        d = float(S.delta(tt))
        worst_d = max(worst_d, abs(xi / d - 1))
        film = x < 0.9 * xi
        eT = np.sqrt(np.mean((T[film] - S.T(x[film], tt)) ** 2)) / (md["Tw"] - md["Tsat"])
        worst_T = max(worst_T, eT)
        # heat the similarity solution conducts to the interface since t0 (all of it evaporates there),
        # against what the solver did with it: evaporation at either face of the slab, and warming it
        q_ref = Lv * md["rho_v"] * (d - md["delta0"])
        ev = m0 - liquid_mass(case, t, md)
        warm = liquid_heat(t) - h0
        worst_e = max(worst_e, abs((Lv * ev + warm) / q_ref - 1))
        rows.append((float(t), xi, d, ev, q_ref / Lv, warm / q_ref, T))
    t, xi, d, ev, mr, fw = rows[-1][:6]
    ok = worst_d < 0.02 and worst_T < 0.02 and worst_e < 0.015
    msg = (f"interface within {worst_d:.2%} of 2 lam sqrt(alpha_v t) (film {md['delta0'] * 1e6:.1f} -> "
           f"{d * 1e6:.1f} um, {xi / d - 1:+.2%} at the end); heat to the slab accounted for within {worst_e:.2%} "
           f"(evaporated {ev:.4e} kg/m2, {ev / mr - 1:+.2%} on the reference, and {fw:.2%} warmed the slab "
           f"by {fw * Lv * mr / (cl * md.get('rho_l', 7030.0) * md['H']) * 1e3:.0f} mK); film T L2 error "
           f"{worst_T:.2%} of (Tw - Tsat)")
    if plt:
        fig, ax = plt.subplots(1, 2, figsize=(9, 3))
        tq = np.linspace(md["t0"], md["t0"] + rows[-1][0], 100)
        ax[0].plot((tq - md["t0"]) * 1e6, S.delta(tq) * 1e6, "k-", label="2 lam sqrt(alpha t)")
        ax[0].plot([r[0] * 1e6 for r in rows], [r[1] * 1e6 for r in rows], "ro", label="alpha = 0.5")
        ax[0].set_xlabel("t - t0 [us]"); ax[0].set_ylabel("film thickness [um]"); ax[0].legend(fontsize=7)
        for r in rows[:: max(1, len(rows) // 4)] + [rows[-1]]:
            tt = md["t0"] + r[0]
            ax[1].plot(x * 1e6, r[6], ".", ms=3)
            xx = np.linspace(0, r[2], 100)
            ax[1].plot(xx * 1e6, S.T(xx, tt), "k-", lw=0.8)
        ax[1].set_xlabel("x [um]"); ax[1].set_ylabel("T [K]")
        savefig(fig, case, "vv_stefanBoiling.png")
    return ("PASS" if ok else "FAIL"), msg


def check_stefan_condensation(case):
    from stefan_lv import KineticCondensation
    md = load(case)
    tl = times(case)
    x = md["x"]
    K = KineticCondensation(k=md["k_l"], rho=md["rho_l"], c=md["c_l"], T0=md["T0"], p_inf=md["p_inf"],
                            cp_v=md["cp_v"], Tsat=md["Tsat"], sigma=md["sigma_c"],
                            t_end=float(tl[-1]) * 1.0001, n=3000)
    m0 = liquid_mass(case, tl[0], md)
    dT = md["Tsat"] - md["T0"]
    rows, worst_m, worst_T = [], 0.0, 0.0
    for t in tl[1:]:
        tt = float(t)
        al, T = F(case, t, "alpha.metal1", md), F(case, t, "T", md)
        # surface: where alpha falls through 0.5 (liquid below, vapour above)
        i = int(np.argmax(al < 0.5))
        xs = float(x[i - 1] + (al[i - 1] - 0.5) / (al[i - 1] - al[i]) * md["dx"])
        mc = liquid_mass(case, t, md) - m0
        mr = float(K.mass(tt))
        if tt >= 2e-6:
            worst_m = max(worst_m, abs(mc / mr - 1))
        liq = (al > 0.99) & (x > xs - 10e-6)
        Tr = K.T(xs - x[liq], tt)
        worst_T = max(worst_T, np.sqrt(np.mean((T[liq] - Tr) ** 2)) / dT)
        rows.append((tt, mc, mr, xs, T, liq, Tr))
    mc, mr = rows[-1][1], rows[-1][2]
    ok = worst_m < 0.03 and worst_T < 0.02
    tt = rows[-1][0]
    msg = (f"condensed mass {mc:.4e} vs {mr:.4e} kg/m2 at {tt * 1e6:.0f} us ({mc / mr - 1:+.2%}; worst after 2 us "
           f"{worst_m:.2%}); rate there {np.interp(tt, K.t, K.j):.2f} kg/m2/s, surface {K.surface(tt):.1f} K "
           f"(from {md['T0']:.0f}); liquid T L2 error {worst_T:.2%} of (Tsat - T0) in the top 10 um")
    if plt:
        fig, ax = plt.subplots(1, 2, figsize=(9, 3))
        ax[0].plot(K.t * 1e6, K.m * 1e3, "k-", label="Hertz-Knudsen + conduction")
        ax[0].plot([r[0] * 1e6 for r in rows], [r[1] * 1e3 for r in rows], "ro", ms=3, label="solver")
        ax[0].set_xlabel("t [us]"); ax[0].set_ylabel("condensed mass [g/m2]"); ax[0].legend(fontsize=7)
        for r in rows[:: max(1, len(rows) // 4)] + [rows[-1]]:
            ax[1].plot((r[3] - x[r[5]]) * 1e6, r[4][r[5]], ".", ms=3)
            ax[1].plot((r[3] - x[r[5]]) * 1e6, r[6], "k-", lw=0.8)
        ax[1].set_xlabel("depth below the surface [um]"); ax[1].set_ylabel("T [K]")
        savefig(fig, case, "vv_stefanCondensation.png")
    return ("PASS" if ok else "FAIL"), msg


def check_thermocapillary(case):
    """19: steady thermocapillary Couette layer. The liquid shear stress must equal the Marangoni stress
    dsigma/dT*G (the tangential stress balance, whatever the interface smearing), the argon must ride at
    the surface speed, and that speed is tau*(h - dx/2)/mu_l (the force sits in the top liquid cell)"""
    md = load(case)
    Nx, Ny, dx = md["Nx"], md["Ny"], md["dx"]
    tl = times(case)
    tau, mu, h = md["tau"], md["mu_l"], md["h"]

    def profile(t):
        U = read_field(case, t, "U", Nx * Ny)
        return U[:, 1].reshape(Ny, Nx).mean(0), np.abs(U[:, 0]).max()

    uy, uxmax = profile(tl[-1])
    uy0, _ = profile(tl[-2])
    xc = (np.arange(Nx) + 0.5) * dx
    nl = int(round(h / dx))
    core = slice(3, nl - 3)
    tau_m = mu * np.polyfit(xc[core], uy[core], 1)[0]
    us_ref = tau * (h - 0.5 * dx) / mu
    u_gas = float(np.mean(uy[nl + 2:]))
    steady = abs(uy[nl - 1] / uy0[nl - 1] - 1)
    e_tau, e_us = tau_m / tau - 1, u_gas / us_ref - 1
    ok = abs(e_tau) < 0.01 and abs(e_us) < 0.01 and steady < 0.005 and uxmax < 1e-3 * abs(us_ref)
    msg = (f"liquid shear stress {tau_m:.2f} Pa vs dsigma/dT*G = {tau:.2f} Pa ({e_tau:+.2%}); argon at "
           f"{u_gas:.4f} m/s vs tau*(h - dx/2)/mu = {us_ref:.4f} ({e_us:+.2%}; tau*h/mu = {tau * h / mu:.4f}); "
           f"change over the last write {steady:.2%}, max|u_n| {uxmax:.1e} m/s")
    if plt:
        fig, ax = plt.subplots(figsize=(4.5, 3))
        ax.plot(uy, xc * 1e6, "ro", ms=3, label="solver")
        xs = np.linspace(0, h, 50)
        ax.plot(tau * xs / mu, xs * 1e6, "k-", lw=0.8, label="Couette, tau = dsigma/dT G")
        ax.axhline(h * 1e6, color="0.6", lw=0.6)
        ax.set_xlabel("u along the surface [m/s]"); ax.set_ylabel("height [um]"); ax.legend(fontsize=7)
        savefig(fig, case, "vv_thermocapillary.png")
    return ("PASS" if ok else "FAIL"), msg


CHECKS = [
    ("01_acousticPulse", "sound speed and adiabatic T-p coupling (chi = 1/gamma closure)", check_acoustic),
    ("02_thermalBox", "energy and closure, one gas: sealed box keeps p0, T -> harmonic mean", check_thermal_box),
    ("03_thermalBoxTwoGas", "energy and closure, two gases (gamma 5/3 and 7/5): exact end state", check_two_gas_box),
    ("04_sod", "shock tube: p*, u*, shock speed (Riemann solution)", check_sod),
    ("05_contactArgon", "advected hot/cold contact, equal gamma: U and p stay uniform", lambda c: check_contact(c, False)),
    ("06_contactAir", "advected contact, gamma 5/3 into 7/5: mixing contraction", lambda c: check_contact(c, True)),
    ("07_evapBoxIsothermal", "sealed evaporation, isothermal walls: p -> Psat(T0)", check_evap_iso),
    ("08_evapBoxAdiabatic", "sealed evaporation, adiabatic: latent heat + Psat end state", check_evap_adiabatic),
    ("09_evapSteady", "steady evaporation: flux, sharp film, recoil jump", check_evap_steady),
    ("10_stefanMelting", "melting front (three-zone similarity solution)", check_stefan),
    ("17_stefanBoiling", "film boiling: film growth and T (Stefan similarity solution), heat budget",
     check_stefan_boiling),
    ("18_stefanCondensation", "condensation onto subcooled liquid: Hertz-Knudsen + conduction (Duhamel)",
     check_stefan_condensation),
    ("19_thermocapillaryLayer", "Marangoni stress: tangential stress balance of a thermocapillary Couette layer",
     check_thermocapillary),
]


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(HERE), "cases")
    rows = []
    for name, what, fn in CHECKS:
        case = os.path.join(root, name)
        if not os.path.isdir(case) or len(times(case)) < 2:
            rows.append((name, what, "NOT RUN", ""))
            continue
        try:
            st, msg = fn(case)
        except Exception as e:                    # report, keep going
            st, msg = "ERROR", f"{type(e).__name__}: {e}"
        rows.append((name, what, st, msg))
        print(f"{st:7s} {name:22s} {msg}")
    with open(os.path.join(root, "vv_results.md"), "w") as f:
        f.write("| case | tests | result | numbers |\n|---|---|---|---|\n")
        for r in rows:
            f.write("| " + " | ".join(r) + " |\n")
    return 0 if all(r[2] in ("PASS", "KNOWN") for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
