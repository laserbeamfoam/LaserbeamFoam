"""ASCII polyMesh reader and cell geometry (numpy only), for checking 3D AMR runs.

Each written time of an AMR run carries its own polyMesh. A decomposed run is read processor by
processor and the cells are concatenated, so nothing needs reconstructing: every check is a sum
or a mean over cells.
"""
import os
import re

import numpy as np

from foamio import read_field, times as root_times

_int = re.compile(r"-?\d+")


def _body(path):
    """Text after the FoamFile header, and the header itself."""
    with open(path) as f:
        s = f.read()
    m = re.search(r"FoamFile\s*\{.*?\}", s, re.S)
    head = m.group(0) if m else ""
    return s[m.end():] if m else s, head


def _list_start(s):
    """Count and position just after the opening '(' of the top-level list."""
    s = re.sub(r"//[^\n]*", "", s)
    m = re.search(r"(\d+)\s*\(", s)
    return int(m.group(1)), s, m.end()


def read_labels(path):
    s, _ = _body(path)
    n, s, i = _list_start(s)
    j = s.index(")", i)
    v = np.array(s[i:j].split(), dtype=np.int64)
    assert v.size == n, path
    return v


def read_points(path):
    s, _ = _body(path)
    n, s, i = _list_start(s)
    j = s.rindex(")")
    v = np.array(s[i:j].replace("(", " ").replace(")", " ").split(), dtype=float)
    assert v.size == 3 * n, path
    return v.reshape(n, 3)


def read_faces(path):
    """Faces as (sizes, flat vertex list); faceList ('4(a b c d)') or faceCompactList."""
    s, head = _body(path)
    if "faceCompactList" in head:
        n, s, i = _list_start(s)
        j = s.index(")", i)
        offs = np.array(s[i:j].split(), dtype=np.int64)
        assert offs.size == n
        m, s2, k = _list_start(s[j + 1:])
        flat = np.array(s2[k:s2.index(")", k)].split(), dtype=np.int64)
        return np.diff(offs), flat
    n, s, i = _list_start(s)
    rows = re.findall(r"(\d+)\(([^)]*)\)", s[i:])
    assert len(rows) == n, path
    sizes = np.array([int(r[0]) for r in rows], dtype=np.int64)
    flat = np.array(" ".join(r[1] for r in rows).split(), dtype=np.int64)
    assert flat.size == sizes.sum(), path
    return sizes, flat


def geometry(points, sizes, flat, owner, neighbour):
    """Cell volumes and centroids by OpenFOAM's own face-triangle / cell-pyramid decomposition."""
    nf = sizes.size
    ncells = int(max(owner.max(), neighbour.max() if neighbour.size else 0)) + 1
    starts = np.concatenate(([0], np.cumsum(sizes)[:-1]))
    fid = np.repeat(np.arange(nf), sizes)

    def fsum(w, ids, n):
        return np.stack([np.bincount(ids, w[:, k], n) for k in range(3)], axis=1)

    pts = points[flat]
    pavg = fsum(pts, fid, nf) / sizes[:, None]
    nxt = np.arange(flat.size) + 1
    nxt[starts + sizes - 1] = starts
    p1, p2, pc = pts, points[flat[nxt]], pavg[fid]
    a = 0.5 * np.cross(p2 - p1, pc - p1)
    Sf = fsum(a, fid, nf)
    nhat = Sf / np.linalg.norm(Sf, axis=1)[:, None]
    am = np.einsum("ij,ij->i", a, nhat[fid])
    Cf = fsum(am[:, None] * (p1 + p2 + pc) / 3.0, fid, nf) / np.bincount(fid, am, nf)[:, None]

    ni = neighbour.size
    nfc = np.bincount(owner, minlength=ncells) + np.bincount(neighbour, minlength=ncells)
    cEst = (fsum(Cf, owner, ncells) + fsum(Cf[:ni], neighbour, ncells)) / nfc[:, None]
    v3o = np.einsum("ij,ij->i", Sf, Cf - cEst[owner])
    v3n = np.einsum("ij,ij->i", Sf[:ni], cEst[neighbour] - Cf[:ni])
    pco = 0.75 * Cf + 0.25 * cEst[owner]
    pcn = 0.75 * Cf[:ni] + 0.25 * cEst[neighbour]
    V3 = np.bincount(owner, v3o, ncells) + np.bincount(neighbour, v3n, ncells)
    C = (fsum(v3o[:, None] * pco, owner, ncells) + fsum(v3n[:, None] * pcn, neighbour, ncells)) / V3[:, None]
    return V3 / 3.0, C


def roots(case):
    """processor*/ directories when the run was decomposed (and not reconstructed), else the case."""
    procs = sorted((d for d in os.listdir(case) if re.fullmatch(r"processor\d+", d)),
                   key=lambda d: int(d[9:]))
    procs = [os.path.join(case, d) for d in procs]
    if procs and len(root_times(procs[0])) > len(root_times(case)):
        return procs
    return [case]


def times(case):
    return root_times(roots(case)[0])


def _mesh_dir(root, time):
    ts = root_times(root)
    earlier = [t for t in ts if float(t) <= float(time)]
    for t in reversed(earlier):
        d = os.path.join(root, t, "polyMesh")
        if os.path.isfile(os.path.join(d, "faces")) or os.path.isfile(os.path.join(d, "faces.gz")):
            return d
    return os.path.join(root, "constant", "polyMesh")


_cache = {}


def mesh(root, time):
    """(V, C) of one root at one time; meshes shared between times are read once."""
    d = _mesh_dir(root, time)
    if d not in _cache:
        pts = read_points(os.path.join(d, "points"))
        sizes, flat = read_faces(os.path.join(d, "faces"))
        own = read_labels(os.path.join(d, "owner"))
        nei = read_labels(os.path.join(d, "neighbour"))
        if len(_cache) > 64:
            _cache.clear()
        _cache[d] = geometry(pts, sizes, flat, own, nei)
    return _cache[d]


class Snapshot:
    """Every root's cells at one time, concatenated: V, C and fields on demand."""

    def __init__(self, case, time):
        self.case, self.time = case, time
        self.parts = []
        for r in roots(case):
            V, C = mesh(r, time)
            self.parts.append((r, V, C))
        self.V = np.concatenate([p[1] for p in self.parts])
        self.C = np.concatenate([p[2] for p in self.parts])
        self.n = self.V.size
        self._f = {}

    def __getitem__(self, name):
        if name not in self._f:
            self._f[name] = np.concatenate([read_field(r, self.time, name, n=V.size)
                                            for r, V, _ in self.parts])
        return self._f[name]

    def has(self, name):
        return all(os.path.isfile(os.path.join(r, self.time, name)) for r, _, _ in self.parts)
