"""Exact Riemann solver for the 1D Euler equations, ideal gas (Toro, ch. 4)."""
import math

import numpy as np


def _f(p, rho, pk, ck, g):
    if p > pk:                                   # shock
        A = 2.0 / ((g + 1) * rho)
        B = (g - 1) / (g + 1) * pk
        return (p - pk) * math.sqrt(A / (p + B))
    return 2 * ck / (g - 1) * ((p / pk) ** ((g - 1) / (2 * g)) - 1)   # rarefaction


def star(rhoL, uL, pL, rhoR, uR, pR, g):
    cL, cR = math.sqrt(g * pL / rhoL), math.sqrt(g * pR / rhoR)
    lo, hi = 1e-12 * min(pL, pR), 100 * max(pL, pR)
    F = lambda p: _f(p, rhoL, pL, cL, g) + _f(p, rhoR, pR, cR, g) + uR - uL
    for _ in range(200):                         # bisection: robust, F is monotone in p
        mid = 0.5 * (lo + hi)
        if F(mid) > 0:
            hi = mid
        else:
            lo = mid
    ps = 0.5 * (lo + hi)
    us = 0.5 * (uL + uR) + 0.5 * (_f(ps, rhoR, pR, cR, g) - _f(ps, rhoL, pL, cL, g))
    return ps, us


def sample(x, t, x0, rhoL, uL, pL, rhoR, uR, pR, g):
    """rho, u, p at positions x (array) and time t > 0."""
    ps, us = star(rhoL, uL, pL, rhoR, uR, pR, g)
    cL, cR = math.sqrt(g * pL / rhoL), math.sqrt(g * pR / rhoR)
    out = np.zeros((3, len(x)))
    for i, xi in enumerate(x):
        S = (xi - x0) / t
        if S <= us:                                         # left of contact
            if ps > pL:                                     # left shock
                SL = uL - cL * math.sqrt((g + 1) / (2 * g) * ps / pL + (g - 1) / (2 * g))
                if S <= SL:
                    out[:, i] = rhoL, uL, pL
                else:
                    r = rhoL * (ps / pL + (g - 1) / (g + 1)) / ((g - 1) / (g + 1) * ps / pL + 1)
                    out[:, i] = r, us, ps
            else:                                           # left rarefaction
                SHL = uL - cL
                cs = cL * (ps / pL) ** ((g - 1) / (2 * g))
                STL = us - cs
                if S <= SHL:
                    out[:, i] = rhoL, uL, pL
                elif S >= STL:
                    out[:, i] = rhoL * (ps / pL) ** (1 / g), us, ps
                else:
                    c = 2 / (g + 1) * (cL + (g - 1) / 2 * (uL - S))
                    u = 2 / (g + 1) * (cL + (g - 1) / 2 * uL + S)
                    r = rhoL * (c / cL) ** (2 / (g - 1))
                    out[:, i] = r, u, pL * (c / cL) ** (2 * g / (g - 1))
        else:                                               # right of contact
            if ps > pR:                                     # right shock
                SR = uR + cR * math.sqrt((g + 1) / (2 * g) * ps / pR + (g - 1) / (2 * g))
                if S >= SR:
                    out[:, i] = rhoR, uR, pR
                else:
                    r = rhoR * (ps / pR + (g - 1) / (g + 1)) / ((g - 1) / (g + 1) * ps / pR + 1)
                    out[:, i] = r, us, ps
            else:
                SHR = uR + cR
                cs = cR * (ps / pR) ** ((g - 1) / (2 * g))
                STR = us + cs
                if S >= SHR:
                    out[:, i] = rhoR, uR, pR
                elif S <= STR:
                    out[:, i] = rhoR * (ps / pR) ** (1 / g), us, ps
                else:
                    c = 2 / (g + 1) * (cR - (g - 1) / 2 * (uR - S))
                    u = 2 / (g + 1) * (-cR + (g - 1) / 2 * uR + S)
                    r = rhoR * (c / cR) ** (2 / (g - 1))
                    out[:, i] = r, u, pR * (c / cR) ** (2 * g / (g - 1))
    return out


def shock_speed(rhoR, uR, pR, ps, g):
    cR = math.sqrt(g * pR / rhoR)
    return uR + cR * math.sqrt((g + 1) / (2 * g) * ps / pR + (g - 1) / (2 * g))
