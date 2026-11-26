"""Minimal readers for ASCII OpenFOAM fields and solver logs (numpy only)."""
import os
import re

import numpy as np

_num = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"


def times(case):
    out = []
    for d in os.listdir(case):
        try:
            t = float(d)
        except ValueError:
            continue
        if os.path.isdir(os.path.join(case, d)):
            out.append((t, d))
    return [d for t, d in sorted(out)]


def read_field(case, time, name, n=None):
    """internalField of a vol field as a numpy array ((n,) scalars or (n,3) vectors)."""
    path = os.path.join(case, time, name)
    with open(path) as f:
        s = f.read()
    m = re.search(r"internalField\s+(uniform|nonuniform)", s)
    if not m:
        raise ValueError(f"no internalField in {path}")
    rest = s[m.end():]
    if m.group(1) == "uniform":
        v = rest.strip()
        if v.startswith("("):
            vec = np.array([float(x) for x in re.findall(_num, v[:v.index(")")])])
            return np.tile(vec, (n or 1, 1))
        val = float(re.match(_num, v).group(0))
        return np.full(n or 1, val)
    m2 = re.search(r"List<(scalar|vector)>\s*(\d+)\s*\(", rest)
    kind, cnt = m2.group(1), int(m2.group(2))
    body = rest[m2.end():]
    if kind == "scalar":
        vals = np.array(re.findall(_num, body[:body.index(")")]), dtype=float)[:cnt]
        return vals
    end = body.index("\n)")
    vals = np.array(re.findall(_num, body[:end]), dtype=float)[:3 * cnt]
    return vals.reshape(cnt, 3)


def log_series(logfile, pattern, groups=1):
    """All matches of a regex in the log, as float arrays (one column per group)."""
    rx = re.compile(pattern)
    rows = []
    with open(logfile, errors="replace") as f:
        for line in f:
            m = rx.search(line)
            if m:
                rows.append([float(m.group(i + 1)) for i in range(groups)])
    return np.array(rows)
