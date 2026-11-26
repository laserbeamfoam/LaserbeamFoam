#!/usr/bin/env python3
"""References for the liquid-vapour Stefan cases (17, 18).

17  Film boiling (one-phase Stefan problem in the vapour, Welch & Wilson 2000): a wall at Tw > Tsat
    under a vapour film, liquid at Tsat beyond it. Conduction in the film feeds the evaporation:
        delta(t) = 2 lam sqrt(alpha_v t),   lam exp(lam^2) erf(lam) = St/sqrt(pi),   St = cp_v (Tw - Tsat)/L
        T(x, t)  = Tw - (Tw - Tsat) erf(x / (2 sqrt(alpha_v t))) / erf(lam)
    The interface sits at Tsat: Hertz-Knudsen and the inertia of the liquid it pushes shift it by
    < 0.3 K here (see gen_cases.py), which the tolerances cover.

18  Condensation onto subcooled liquid with Hertz-Knudsen kinetics at the surface. Saturated vapour at
    p_inf, Tsat over liquid at T0 < Tsat. For liquid Fe the kinetics, not conduction, limit the rate
    for t << (k rho c)/(L kappa dPsat/dT)^2 ~ 1 ms, so the classical Tsat-surface solution does not
    apply. The liquid is a half-space (it moves by < 0.1 um), the surface flux is
        q(t) = (L + cp_v (Tsat - T_s)) j(T_s),   j(T) = sigma sqrt(M/(2 pi R T)) (p_inf - Psat(T)),
    and Duhamel gives the surface temperature
        T_s(t) = T0 + 1/sqrt(pi k rho c) int_0^t q(s)/sqrt(t - s) ds,
    a nonlinear Abel-Volterra equation, solved here by product integration (q linear on each step,
    exact weights, Newton on the newest node). Interior temperatures follow from the same q(t).
"""
import math

import numpy as np
from scipy.optimize import brentq
from scipy.special import erf

RU = 8.314            # the solver's constant in Clausius-Clapeyron and Hertz-Knudsen
ML, LV, TB, P0 = 0.05585, 6.09e6, 3134.0, 1e5


def psat(T):
    return P0 * np.exp(ML * LV / (TB * RU) * (1 - TB / np.asarray(T, dtype=float)))


def dpsat(T):
    return psat(T) * ML * LV / (RU * np.asarray(T, dtype=float) ** 2)


def j_hk(T, p, sigma=1.0):
    """Hertz-Knudsen mass flux [kg/m2/s], positive = condensation onto liquid at T from vapour at p"""
    return sigma * np.sqrt(ML / (2 * math.pi * RU * T)) * (p - psat(T))


class FilmBoiling:
    def __init__(self, k, rho, cp, Tw, Tsat=TB, L=LV):
        self.alpha = k / (rho * cp)
        self.Tw, self.Tsat, self.rho, self.k, self.L = Tw, Tsat, rho, k, L
        self.St = cp * (Tw - Tsat) / L
        self.lam = brentq(lambda a: a * math.exp(a * a) * math.erf(a) - self.St / math.sqrt(math.pi), 1e-9, 5.0)

    def delta(self, t):
        return 2 * self.lam * np.sqrt(self.alpha * np.asarray(t, dtype=float))

    def T(self, x, t):
        x = np.asarray(x, dtype=float)
        return self.Tw - (self.Tw - self.Tsat) * erf(x / (2 * np.sqrt(self.alpha * t))) / math.erf(self.lam)

    def flux(self, t):
        """evaporation mass flux [kg/m2/s] = rho_v d(delta)/dt"""
        return self.rho * self.lam * np.sqrt(self.alpha / np.asarray(t, dtype=float))


class KineticCondensation:
    def __init__(self, k, rho, c, T0, p_inf, cp_v, Tsat=TB, L=LV, sigma=1.0, t_end=1e-4, n=4000):
        self.k, self.rho, self.c, self.T0 = k, rho, c, T0
        self.alpha = k / (rho * c)
        self.e = math.sqrt(math.pi * k * rho * c)
        self.p, self.cpv, self.Tsat, self.L, self.sig = p_inf, cp_v, Tsat, L, sigma
        self._solve(t_end, n)

    def q_of(self, Ts):
        return (self.L + self.cpv * (self.Tsat - Ts)) * j_hk(Ts, self.p, self.sig)

    def _solve(self, t_end, n):
        # graded steps (fine where q changes fastest), product integration with q linear per step
        t = t_end * (np.linspace(0.0, 1.0, n + 1) ** 2)
        Ts = np.empty(n + 1)
        q = np.empty(n + 1)
        Ts[0] = self.T0
        q[0] = self.q_of(self.T0)
        for m in range(1, n + 1):
            tm = t[m]
            # int_{t_i}^{t_i+1} (linear q)/sqrt(tm - s) ds: weights on q_i and q_i+1
            a, b = tm - t[:m], tm - t[1:m + 1]
            sa, sb = np.sqrt(a), np.sqrt(b)
            h = t[1:m + 1] - t[:m]
            I0 = 2 * (sa - sb)                                   # int 1/sqrt
            I1 = (2.0 / 3.0) * (a * sa - b * sb)                 # int (tm - s)/sqrt
            # q(s) = q_i + (q_i+1 - q_i)(s - t_i)/h,  s - t_i = a - (tm - s)
            w_i = I0 - (a * I0 - I1) / h
            w_ip1 = (a * I0 - I1) / h
            known = np.dot(w_i, q[:m]) + np.dot(w_ip1[:-1], q[1:m])
            wl = w_ip1[-1]

            def f(T):
                return T - self.T0 - (known + wl * self.q_of(T)) / self.e

            Ts[m] = brentq(f, self.T0 - 50.0, self.Tsat + 50.0)
            q[m] = self.q_of(Ts[m])
        self.t, self.Ts, self.q = t, Ts, q
        self.j = q / (self.L + self.cpv * (self.Tsat - Ts))
        # condensed mass per unit area [kg/m2], trapezoidal in the graded steps
        self.m = np.concatenate(([0.0], np.cumsum(0.5 * (self.j[1:] + self.j[:-1]) * np.diff(t))))

    def surface(self, tq):
        return np.interp(tq, self.t, self.Ts)

    def mass(self, tq):
        return np.interp(tq, self.t, self.m)

    def T(self, z, tq):
        """temperature at depth z below the surface at time tq (half-space, same q(t))"""
        z = np.atleast_1d(np.asarray(z, dtype=float))
        out = np.full(z.shape, self.T0)
        s = self.t[self.t < tq]
        if len(s) < 2:
            return out
        qs = np.interp(s, self.t, self.q)
        # q(s) G(z, tq - s), G = exp(-z^2/(4 alpha tau))/sqrt(tau); midpoint on the graded steps
        sm = 0.5 * (s[1:] + s[:-1])
        qm = 0.5 * (qs[1:] + qs[:-1])
        ds = np.diff(s)
        tau = tq - sm
        for i, zi in enumerate(z):
            out[i] += np.sum(qm * np.exp(-zi * zi / (4 * self.alpha * tau)) / np.sqrt(tau) * ds) / self.e
        # the last partial step to tq (q ~ q(s[-1])), integrable 1/sqrt singularity at z = 0
        tail = tq - s[-1]
        if tail > 0:
            out += qs[-1] * 2 * np.sqrt(tail) * np.exp(-z * z / (4 * self.alpha * max(tail, 1e-30))) / self.e
        return out


if __name__ == "__main__":
    fb = FilmBoiling(k=1e-5 * 2.5 * 8.314 / 0.05585 / 0.0067, rho=1e5 * 55.85 / (8314.47 * TB),
                     cp=2.5 * 8.314 / 0.05585, Tw=TB + 50.0)
    print(f"film boiling: St {fb.St:.4e}, lambda {fb.lam:.6f}, alpha {fb.alpha:.4e} m2/s, "
          f"delta(100 us) {fb.delta(1e-4) * 1e6:.2f} um, delta(400 us) {fb.delta(4e-4) * 1e6:.2f} um, "
          f"j(100 us) {fb.flux(1e-4):.4f} kg/m2/s")
    kc = KineticCondensation(k=30.0, rho=7030.0, c=800.0, T0=TB - 100.0, p_inf=1e5,
                             cp_v=2.5 * 8.314 / 0.05585, t_end=2e-5)
    for tq in (1e-7, 1e-6, 5e-6, 1e-5, 2e-5):
        print(f"condensation t {tq * 1e6:5.1f} us: T_s {kc.surface(tq):8.2f} K, j {np.interp(tq, kc.t, kc.j):6.3f} "
              f"kg/m2/s, mass {kc.mass(tq):.4e} kg/m2, T at 0.25 um {kc.T(0.25e-6, tq)[0]:.2f}")
