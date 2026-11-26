#!/usr/bin/env python3
"""
stefan_reference.py -- exact similarity solutions of the one-dimensional
Stefan problem, for verifying the melting/solidification machinery of
laserbeamFoam and compressibleLaserbeamFoam (the epsilon1 correction loop,
both 'sourceBased' and 'apparentCapacity').

Geometry: a semi-infinite bar x >= 0, initially at Ti, whose end x = 0 is held
at Tw for t > 0. Constant rho, k, c in every phase (the same values), latent
heat L.

  pure   (Tsolidus == Tliquidus == Tm): the classical Neumann solution,
         one front at x = 2*lambda*sqrt(t) with a latent-heat jump condition.
  mushy  (Tsolidus < Tliquidus, linear liquid fraction in between): three
         zones -- wall-side phase, mushy zone, far-side phase -- each with
         constant diffusivity, the mushy one with the apparent heat capacity
         c + L/(Tliquidus - Tsolidus). Two fronts x1 = 2*lambda1*sqrt(t)
         (solidus or liquidus, wall side) and x2 = 2*lambda2*sqrt(t);
         T and k*dT/dx are continuous at both (the latent heat is distributed
         over the mushy zone, so there is no flux jump). This is EXACTLY the
         problem the solver's enthalpy-porosity loop solves, so the comparison
         should converge with the mesh to round-off, not to "model error".

Melting   : Tw > Tliquidus, Ti < Tsolidus   (wall-side zone is liquid)
Solidifying: Tw < Tsolidus, Ti > Tliquidus  (wall-side zone is solid)

Usage
  Reference only:
    stefan_reference.py --rho 4430 --k 20 --c 750 --L 2.9e5 --Ts 1878 --Tl 1933 \
        --Tw 2100 --Ti 1500 --times 0.001 0.01 0.1 --out stefan_ref.csv
  Compare with a case (1D along --axis, cell centres from writeCellCentres):
    stefan_reference.py ... --case . --axis x --tol-front 0.05 --tol-T 0.02
    The simulated front is where epsilon1 crosses 0.5 (linear interpolation
    along the bar); the reference front is where the mushy-zone temperature
    equals (Ts+Tl)/2. The temperature error is the L2 norm over the bar
    relative to |Tw - Ti|.

Exit code: 0 pass/reference-only, 1 fail, 2 error.
"""

import argparse
import json
import math
import os
import sys

import numpy as np

try:
    from scipy.special import erf, erfc, erfinv
    from scipy.optimize import fsolve
except ImportError:  # pragma: no cover
    print("stefan_reference.py needs scipy (pip install scipy)", file=sys.stderr)
    sys.exit(2)


def _trapz(y, x):
    f = getattr(np, 'trapezoid', None) or getattr(np, 'trapz')
    return f(y, x)


class Stefan:
    def __init__(self, rho, k, c, L, Ts, Tl, Tw, Ti):
        self.rho, self.k, self.c, self.L = rho, k, c, L
        self.Ts, self.Tl, self.Tw, self.Ti = Ts, Tl, Tw, Ti
        self.alpha = k / (rho * c)
        self.melting = Tw > Ti
        if self.melting and not (Tw > Tl and Ti < Ts):
            raise ValueError("melting needs Tw > Tliquidus and Ti < Tsolidus")
        if not self.melting and not (Tw < Ts and Ti > Tl):
            raise ValueError("solidifying needs Tw < Tsolidus and Ti > Tliquidus")
        # a mushy range below 1e-4 of the driving temperature difference is
        # numerically indistinguishable from a pure melting point (and the
        # three-zone solve is ill-conditioned there): use Neumann at (Ts+Tl)/2
        self.pure = (Tl - Ts) <= 1e-4 * abs(Tw - Ti)
        if self.pure:
            Tm = 0.5 * (Ts + Tl)
            self.Ts = self.Tl = Tm
            self._solve_pure()
        else:
            self._solve_mushy()

    # ---- pure (Neumann) ---------------------------------------------------
    def _solve_pure(self):
        Tm = self.Ts
        a = self.alpha
        St = self.L / self.c

        def F(u):
            u = float(np.asarray(u).ravel()[0])
            if self.melting:
                return ((self.Tw - Tm) / erf(u) - (Tm - self.Ti) / erfc(u)
                        - St * u * math.sqrt(math.pi) * math.exp(u * u))
            # solid on the wall side: the signs flip
            return ((Tm - self.Tw) / erf(u) - (self.Ti - Tm) / erfc(u)
                    - St * u * math.sqrt(math.pi) * math.exp(u * u))
        u = float(fsolve(F, 0.3)[0])
        self.lam = u * math.sqrt(a)          # X(t) = 2*lam*sqrt(t)
        self.u = u

    # ---- mushy (three zones) ----------------------------------------------
    def _solve_mushy(self):
        c, L = self.c, self.L
        aA = self.alpha
        cM = c + L / (self.Tl - self.Ts)
        aM = self.k / (self.rho * cM)
        aB = self.alpha
        self.aA, self.aM, self.aB = aA, aM, aB
        if self.melting:
            T1, T2 = self.Tl, self.Ts    # wall-side front (liquidus), far front (solidus)
        else:
            T1, T2 = self.Ts, self.Tl
        self.T1, self.T2 = T1, T2
        Tw, Ti = self.Tw, self.Ti

        def coeffs(l1, l2):
            e1 = erf(l1 / math.sqrt(aM))
            e2 = erf(l2 / math.sqrt(aM))
            de = e2 - e1
            if abs(de) < 1e-300:
                de = 1e-300 if de >= 0 else -1e-300
            B = (T2 - T1) / de
            A = T1 - B * e1
            return A, B

        def F(x):
            l1, l2 = float(x[0]), float(x[1])
            if l2 <= l1 or l1 <= 0:
                return [1e3, 1e3]
            A, B = coeffs(l1, l2)
            # flux continuity at x1: zone A vs zone M
            gA = (T1 - Tw) / erf(l1 / math.sqrt(aA)) * math.exp(-l1 * l1 / aA) / math.sqrt(aA)
            gM1 = B * math.exp(-l1 * l1 / aM) / math.sqrt(aM)
            # flux continuity at x2: zone M vs zone B
            gM2 = B * math.exp(-l2 * l2 / aM) / math.sqrt(aM)
            gB = -(T2 - Ti) / erfc(l2 / math.sqrt(aB)) * math.exp(-l2 * l2 / aB) / math.sqrt(aB)
            return [gA - gM1, gM2 - gB]

        # initial guess from the pure solution with Tm = (Ts+Tl)/2: lambda1 ~ g,
        # mushy thickness from the pure-front gradient. Unknowns (lambda1, log(delta))
        # keep lambda2 > lambda1; residuals are scaled by the pure-front gradient so
        # the solve stays conditioned for a vanishingly thin mushy range.
        tmp = Stefan(self.rho, self.k, c, L, 0.5 * (self.Ts + self.Tl), 0.5 * (self.Ts + self.Tl), Tw, Ti)
        g = tmp.lam
        u = tmp.u
        grad0 = abs(Tw - Ti) * math.exp(-u * u) / (erf(u) * math.sqrt(math.pi * aA))  # |dT/dx|*sqrt(t) scale
        dlam0 = max(abs(self.Tl - self.Ts) / (2.0 * grad0), 1e-12 * g)

        def G(y):
            l1 = float(y[0])
            l2 = l1 + math.exp(float(y[1]))
            r = F([l1, l2])
            return [r[0] / grad0, r[1] / grad0]

        best = None
        for y0 in ([g - 0.5 * dlam0, math.log(dlam0)], [0.9 * g, math.log(max(dlam0, 0.05 * g))],
                   [0.7 * g, math.log(max(dlam0, 0.3 * g))]):
            sol, info, ier, msg = fsolve(G, y0, full_output=True, xtol=1e-13)
            res = max(abs(v) for v in G(sol))
            if best is None or res < best[0]:
                best = (res, sol, msg)
            if ier == 1 and res < 1e-9:
                break
        res, sol, msg = best
        if res > 1e-6:
            raise RuntimeError(f"mushy Stefan: front solve failed (residual {res:.2e}): {msg}")
        self.lam1 = float(sol[0])
        self.lam2 = self.lam1 + math.exp(float(sol[1]))
        self.A, self.B = coeffs(self.lam1, self.lam2)

    # ---- evaluation ---------------------------------------------------------
    def fronts(self, t):
        """(x_wallside, x_far, x_half) : solidus/liquidus fronts and the eps=0.5 position."""
        s = math.sqrt(max(t, 0.0))
        if self.pure:
            X = 2 * self.lam * s
            return X, X, X
        x1 = 2 * self.lam1 * s
        x2 = 2 * self.lam2 * s
        Tmid = 0.5 * (self.Ts + self.Tl)
        if x2 - x1 < 1e-9 * max(x2, 1e-300) or not np.isfinite(self.B) or self.B == 0:
            return x1, x2, 0.5 * (x1 + x2)
        arg = (Tmid - self.A) / self.B
        xh = 2 * math.sqrt(self.aM * t) * erfinv(min(max(arg, -1 + 1e-15), 1 - 1e-15))
        if not np.isfinite(xh):
            xh = 0.5 * (x1 + x2)
        return x1, x2, xh

    def T(self, x, t):
        x = np.asarray(x, dtype=float)
        out = np.empty_like(x)
        s = math.sqrt(max(t, 1e-300))
        if self.pure:
            Tm = self.Ts
            X = 2 * self.lam * s
            u = self.u
            wall = x <= X
            out[wall] = self.Tw + (Tm - self.Tw) * erf(x[wall] / (2 * math.sqrt(self.alpha * t))) / erf(u)
            out[~wall] = self.Ti + (Tm - self.Ti) * erfc(x[~wall] / (2 * math.sqrt(self.alpha * t))) / erfc(u)
            return out
        x1, x2, _ = self.fronts(t)
        zA = x <= x1
        zB = x >= x2
        zM = ~zA & ~zB
        out[zA] = self.Tw + (self.T1 - self.Tw) * erf(x[zA] / (2 * math.sqrt(self.aA * t))) / erf(self.lam1 / math.sqrt(self.aA))
        out[zM] = self.A + self.B * erf(x[zM] / (2 * math.sqrt(self.aM * t)))
        out[zB] = self.Ti + (self.T2 - self.Ti) * erfc(x[zB] / (2 * math.sqrt(self.aB * t))) / erfc(self.lam2 / math.sqrt(self.aB))
        return out

    def eps(self, x, t):
        T = self.T(x, t)
        if self.pure:
            return np.where(T > self.Ts, 1.0, 0.0)
        return np.clip((T - self.Ts) / (self.Tl - self.Ts), 0.0, 1.0)

    def wall_flux(self, t):
        """Heat flux INTO the bar at x=0, W/m^2 (positive when heating the bar)."""
        s = math.sqrt(t)
        if self.pure:
            dTdx = (self.Ts - self.Tw) / erf(self.u) / math.sqrt(math.pi * self.alpha * t)
        else:
            dTdx = (self.T1 - self.Tw) / erf(self.lam1 / math.sqrt(self.aA)) / math.sqrt(math.pi * self.aA * t)
        return -self.k * dTdx

    def enthalpy(self, t, xmax_factor=12.0, n=200001):
        """Integral of rho*(c*(T-Ti) + L*(eps-eps_i)) dx over the bar [J/m^2]."""
        if t <= 0:
            return 0.0
        xmax = xmax_factor * 2 * math.sqrt(self.alpha * t) + (self.fronts(t)[1] if not self.pure else 0)
        x = np.linspace(0, xmax, n)
        T = self.T(x, t)
        if self.pure:
            e = np.where(T > self.Ts, 1.0, 0.0)
            ei = 1.0 if self.Ti > self.Ts else 0.0
        else:
            e = np.clip((T - self.Ts) / (self.Tl - self.Ts), 0, 1)
            ei = float(np.clip((self.Ti - self.Ts) / (self.Tl - self.Ts), 0, 1))
        return float(_trapz(self.rho * (self.c * (T - self.Ti) + self.L * (e - ei)), x))


def front_from_case(case, tdir, axis, eps_name="epsilon1"):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import foamfield as ff
    C = ff.read_centres(case, tdir)
    x = C[:, "xyz".index(axis)]
    e = ff.read_scalar(case, tdir, eps_name, C.shape[0])
    o = np.argsort(x)
    x, e = x[o], e[o]
    # collapse duplicates across the transverse direction (1D cases are one cell thick)
    xu, idx = np.unique(np.round(x, 15), return_inverse=True)
    eu = np.array([e[idx == i].mean() for i in range(len(xu))])
    return xu, eu, C, o


def crossing(x, f, level):
    """First x where f crosses 'level' (linear interpolation); nan if never."""
    for i in range(len(x) - 1):
        a, b = f[i] - level, f[i + 1] - level
        if a == 0:
            return x[i]
        if a * b < 0:
            return x[i] + (x[i + 1] - x[i]) * a / (a - b)
    return float("nan")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rho", type=float, required=True)
    ap.add_argument("--k", type=float, required=True)
    ap.add_argument("--c", type=float, required=True, help="sensible specific heat (cp or cv, as the solver uses)")
    ap.add_argument("--L", type=float, required=True, help="latent heat of fusion [J/kg]")
    ap.add_argument("--Ts", type=float, required=True)
    ap.add_argument("--Tl", type=float, required=True)
    ap.add_argument("--Tw", type=float, required=True, help="wall temperature at x=0")
    ap.add_argument("--Ti", type=float, required=True, help="initial temperature")
    ap.add_argument("--times", type=float, nargs="*", default=[])
    ap.add_argument("--out", default=None, help="CSV with the fronts at --times (and at the case times)")
    ap.add_argument("--profile-out", default=None, help="CSV with T(x) at --times on a fine grid")
    ap.add_argument("--xmax", type=float, default=None, help="bar length for --profile-out")
    ap.add_argument("--case", default=None, help="compare with this case")
    ap.add_argument("--axis", choices=["x", "y", "z"], default="x")
    ap.add_argument("--eps-name", default="epsilon1")
    ap.add_argument("--T-name", default="T")
    ap.add_argument("--tol-front", type=float, default=0.05, help="relative tolerance on the eps=0.5 front position")
    ap.add_argument("--tol-T", type=float, default=0.02, help="relative L2 tolerance on T (scaled by |Tw-Ti|)")
    ap.add_argument("--json", default=None)
    ap.add_argument("--no-fail", action="store_true")
    args = ap.parse_args()

    S = Stefan(args.rho, args.k, args.c, args.L, args.Ts, args.Tl, args.Tw, args.Ti)
    kind = "pure (Neumann)" if S.pure else "mushy (three-zone)"
    if S.pure:
        print(f"stefan_reference: {kind}, lambda = {S.lam:.6e} m/s^0.5  (X = 2*lambda*sqrt(t))")
    else:
        print(f"stefan_reference: {kind}, lambda1 = {S.lam1:.6e}, lambda2 = {S.lam2:.6e} m/s^0.5")
        # energy self-check: dH/dt must equal the wall flux (both ~ 1/sqrt(t))
        t0 = 1.0
        H = S.enthalpy(t0)
        q = S.wall_flux(t0)
        # H(t) = 2 q(t) t for similarity solutions (q ~ t^-1/2, H ~ t^1/2)
        print(f"  energy self-check at t=1 s: H = {H:.6e} J/m^2, 2*q*t = {2*q*t0:.6e} J/m^2, "
              f"rel. mismatch {abs(H-2*q*t0)/abs(H):.2e}")

    rows = []
    for t in args.times:
        x1, x2, xh = S.fronts(t)
        rows.append((t, x1, x2, xh))

    result = {"kind": kind, "status": "reference-only"}
    failures = []

    if args.case:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import foamfield as ff
        comp = []
        for t, tdir in ff.time_dirs(args.case, skip_zero=True):
            xu, eu, C, o = front_from_case(args.case, tdir, args.axis, args.eps_name)
            xsim = crossing(xu, eu, 0.5)
            x1, x2, xh = S.fronts(t)
            errF = abs(xsim - xh) / max(abs(xh), 1e-30) if xh > 0 else float("nan")
            # temperature profile error
            errT = float("nan")
            pT = os.path.join(args.case, tdir, args.T_name)
            if os.path.exists(pT):
                Tsim = ff.read_scalar(args.case, tdir, args.T_name, C.shape[0])
                xcell = C[:, "xyz".index(args.axis)]
                Tref = S.T(xcell, t)
                errT = float(np.sqrt(np.mean((Tsim - Tref) ** 2)) / abs(args.Tw - args.Ti))
            comp.append((t, xsim, xh, errF, errT))
            rows.append((t, x1, x2, xh))
            okF = (not np.isnan(errF)) and errF <= args.tol_front
            okT = np.isnan(errT) or errT <= args.tol_T
            print(f"  t={t:.6g}: front sim {xsim:.6e} m  ref {xh:.6e} m  rel.err {errF:.3e} {'ok' if okF else 'FAIL'}"
                  f"   T L2 rel.err {errT:.3e} {'ok' if okT else 'FAIL'}")
            if not okF:
                failures.append(f"t={t:g}: front error {errF:.3e} > {args.tol_front}")
            if not okT:
                failures.append(f"t={t:g}: T error {errT:.3e} > {args.tol_T}")
        result["comparison"] = [dict(zip(["t", "xsim", "xref", "frontRelErr", "TRelL2"], c)) for c in comp]
        result["status"] = "FAIL" if failures else "pass"
        result["failures"] = failures

    if args.out:
        rows.sort()
        with open(args.out, "w") as f:
            f.write("# time[s],x_front_wallside[m],x_front_far[m],x_eps0.5[m]\n")
            for r in rows:
                f.write(",".join(f"{v:.9e}" for v in r) + "\n")
        print(f"stefan_reference: wrote {args.out}")

    if args.profile_out and args.times:
        xmax = args.xmax or 8 * 2 * math.sqrt(S.alpha * max(args.times))
        x = np.linspace(0, xmax, 2001)
        with open(args.profile_out, "w") as f:
            f.write("# x[m]," + ",".join(f"T(t={t:g})" for t in args.times) + "\n")
            Ts = [S.T(x, t) for t in args.times]
            for i in range(len(x)):
                f.write(f"{x[i]:.9e}," + ",".join(f"{T[i]:.6f}" for T in Ts) + "\n")
        print(f"stefan_reference: wrote {args.profile_out}")

    if args.json:
        with open(args.json, "w") as f:
            json.dump(result, f, indent=2)

    for fl in failures:
        print(f"  FAILURE: {fl}")
    if args.case:
        print(f"stefan_reference: {'FAIL' if failures else 'PASS'}")
    return 0 if (args.no_fail or not failures) else 1


if __name__ == "__main__":
    sys.exit(main())
