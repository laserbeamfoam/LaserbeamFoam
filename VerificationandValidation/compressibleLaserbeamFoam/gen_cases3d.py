#!/usr/bin/env python3
"""
Generate the 3D verification cases (AMR, no laser) for compressibleLaserbeamFoam.

Materials, numerics and the laser-off set-up are those of gen_cases.py. Each case is one hex
block of cubic cells, refined at run time by dynamicRefineFvMesh.

Initial fields must be set on the refined mesh: AMR copies a cell's values into its children,
so fields set on the base mesh alone leave a mixed band one coarse cell thick where the interface
should be. The case Allrun therefore lets the solver refine the base-mesh fields for
maxRefinement steps, then calls this script with --init. That moves the refined mesh (with its
refinement history, so it can still unrefine) to t = 0 and writes the fields on it, with exact
volume fractions sub-sampled in every cell.

The base-mesh fields are kept in 0.orig. A copy or checkout can drop that folder (a *.orig ignore
rule, as many global gitignores have for merge backups), so the case Allrun rewrites it with --orig
when it is missing.

Usage: python3 gen_cases3d.py [outdir]        write every case (fields on the base mesh)
       python3 gen_cases3d.py --setup CASE    rewrite one case's whole set-up (the top-level Allrun does this
                                              for a case that lost 0.orig or has an older Allrun)
       python3 gen_cases3d.py --orig CASE     rewrite CASE/0.orig (the case Allrun does this if it is missing)
       python3 gen_cases3d.py --init CASE     from the case Allrun, after the refining steps
"""
import json
import math
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MODE = sys.argv[1] if len(sys.argv) > 2 and sys.argv[1] in ("--init", "--orig", "--setup") else None
TARGET = sys.argv[2] if MODE else None
OUT = "cases" if MODE else (sys.argv[1] if len(sys.argv) > 1 else "cases")
sys.argv = sys.argv[:1] + ([OUT] if not MODE else [])   # gen_cases reads its OUT from argv[1]
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "tools"))
import gen_cases as g                         # noqa: E402

SIDES = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")
FACES = {"xmin": "(0 4 7 3)", "xmax": "(1 2 6 5)", "ymin": "(0 1 5 4)",
         "ymax": "(3 7 6 2)", "zmin": "(0 3 2 1)", "zmax": "(4 5 6 7)"}
OPPOSITE = {"xmin": "xmax", "xmax": "xmin", "ymin": "ymax", "ymax": "ymin", "zmin": "zmax", "zmax": "zmin"}


# ----------------------------------------------------------------------------- geometry
def centres(n, L):
    """Cell centres of the blockMesh block, in blockMesh order (x fastest)."""
    nx, ny, nz = n
    dx = [L[i] / n[i] for i in range(3)]
    k, j, i = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx), indexing="ij")
    return np.stack([(i.ravel() + 0.5) * dx[0], (j.ravel() + 0.5) * dx[1], (k.ravel() + 0.5) * dx[2]], 1)


def fraction(C, h, inside, sub=8):
    """Fraction of each cubic cell (centre C, edge h) where inside(x, y, z) holds, sub^3 samples."""
    o = (np.arange(sub) + 0.5) / sub - 0.5
    ox, oy, oz = np.meshgrid(o, o, o, indexing="ij")
    frac = np.zeros(len(C))
    for s in np.stack([ox.ravel(), oy.ravel(), oz.ravel()], 1):
        p = C + s[None, :] * h[:, None]
        frac += inside(p[:, 0], p[:, 1], p[:, 2])
    return frac / sub ** 3


# ----------------------------------------------------------------------------- boundary conditions
def wall(T=None):
    return {"_type": "wall",
            "U": "type fixedValue; value uniform (0 0 0);",
            "p_rgh": "type fixedFluxPressure; value $internalField;",
            "p": "type calculated; value $internalField;",
            "T": f"type fixedValue; value uniform {T!r};" if T else "type zeroGradient;",
            "Laser_boundary": "type fixedValue; value uniform -1;"}


SYMMETRY = {"_type": "symmetryPlane", "*": "type symmetryPlane;"}
CYCLIC = {"_type": "cyclic", "*": "type cyclic;"}


def bc_text(kind, field, uniform_value):
    txt = kind.get(field, kind.get("*", "type zeroGradient;"))
    return txt.replace("$internalField", uniform_value)


# ----------------------------------------------------------------------------- fields
def write_fields(case, spec, f, tdir="0"):
    """tdir/ for the field values f (arrays over the cells, U a uniform vector)."""
    sides, N = spec["sides"], len(f["T"])
    U0 = f.get("U", (0.0, 0.0, 0.0))

    def one(fname, vals, dims, vector=False):
        if vector:
            body_int = uni = "uniform ({!r} {!r} {!r})".format(*(float(u) for u in U0))
        else:
            vals = np.asarray(vals, dtype=float)
            if np.all(vals == vals[0]):
                body_int = f"uniform {float(vals[0])!r}"
            else:
                body_int = "nonuniform List<scalar>\n{}\n(\n{}\n)".format(N, "\n".join(repr(float(v)) for v in vals))
            uni = f"uniform {float(np.mean(vals))!r}"
        bfs = "\n".join(f"    {s:5s} {{ {bc_text(sides[s], fname, uni)} }}" for s in SIDES)
        g.write(case, f"{tdir}/{fname}", "volVectorField" if vector else "volScalarField", f"""dimensions      {dims};

internalField   {body_int};

boundaryField
{{
{bfs}
}}
""")

    one("T", f["T"], "[0 0 0 1 0 0 0]")
    one("p", f["p"], "[1 -1 -2 0 0 0 0]")
    one("p_rgh", f["p"], "[1 -1 -2 0 0 0 0]")
    one("U", None, "[0 1 -1 0 0 0 0]", vector=True)
    for ph in ("metal1", "metal1vapour", spec["gas"]):
        one(f"alpha.{ph}", f[f"alpha.{ph}"], "[0 0 0 0 0 0 0]")
    one("Laser_boundary", np.zeros(N), "[0 0 0 0 0 0 0]")


# ----------------------------------------------------------------------------- case writer
ALLRUN = """#!/bin/sh
# serial by default; NPROCS=8 ./Allrun decomposes after the initialisation (fields are read per
# processor by tools/check3d.py, so nothing needs reconstructing)
cd "${0%/*}" || exit 1
SOLVER=${SOLVER:-compressibleLaserbeamFoam}
GEN=${GEN:-../../gen_cases3d.py}
name=$(basename "$PWD")
[ -f system/controlDict.run ] && mv system/controlDict.run system/controlDict   # an interrupted earlier run
rm -rf 0 [1-9]* 0.[0-9]* processor* postProcessing log.*
# base-mesh fields (gen_cases3d.py --init rewrites 0/ later); rebuilt if a copy left 0.orig behind
if [ ! -d 0.orig ]; then
    python3 "$GEN" --orig . > log.orig 2>&1 || { echo "$name: 0.orig missing and could not be rebuilt (log.orig)"; exit 1; }
fi
cp -r 0.orig 0
blockMesh > log.blockMesh 2>&1 || { echo "$name: blockMesh failed (log.blockMesh)"; exit 1; }

# the solver refines the base-mesh fields for MAXREF steps; gen_cases3d.py --init then moves that
# mesh to t = 0 and writes the fields on it, so the interface starts sharp on the finest cells
cp system/controlDict system/controlDict.run
sed -i -e "s/^stopAt .*/stopAt          nextWrite;/" -e "s/^writeControl .*/writeControl    timeStep;/" \\
       -e "s/^writeInterval .*/writeInterval   MAXREF;/" system/controlDict
$SOLVER > log.prerefine 2>&1
st=$?
mv system/controlDict.run system/controlDict
[ "$st" = 0 ] || { echo "$name: refining steps failed, exit $st (log.prerefine)"; exit "$st"; }
python3 "$GEN" --init . > log.init 2>&1 || { echo "$name: initialisation failed (log.init)"; exit 1; }

NP=${NPROCS:-1}
if [ "$NP" -gt 1 ]; then
    sed -i "s/^numberOfSubdomains.*/numberOfSubdomains $NP;/" system/decomposeParDict
    decomposePar -force > log.decomposePar 2>&1 || { echo "$name: decomposePar failed (log.decomposePar)"; exit 1; }
    mpirun -np "$NP" $SOLVER -parallel > log.solver 2>&1
else
    $SOLVER > log.solver 2>&1
fi
st=$?
echo "$name: solver exit $st"
exit "$st"
"""


def make_case3d(spec):
    name, n, L, gas, sides = spec["name"], spec["n"], spec["L"], spec["gas"], spec["sides"]
    case = os.path.join(OUT, name)
    N = n[0] * n[1] * n[2]

    # --- mesh
    Lx, Ly, Lz = L
    pl = []
    for s in SIDES:
        k = sides[s]
        if k["_type"] == "cyclic":
            pl.append(f"    {s} {{ type cyclic; neighbourPatch {OPPOSITE[s]}; faces ({FACES[s]}); }}")
        else:
            pl.append(f"    {s} {{ type {k['_type']}; faces ({FACES[s]}); }}")
    g.write(case, "system/blockMeshDict", "dictionary", f"""scale   1;

vertices
(
    (0 0 0) ({Lx!r} 0 0) ({Lx!r} {Ly!r} 0) (0 {Ly!r} 0)
    (0 0 {Lz!r}) ({Lx!r} 0 {Lz!r}) ({Lx!r} {Ly!r} {Lz!r}) (0 {Ly!r} {Lz!r})
);

blocks ( hex (0 1 2 3 4 5 6 7) ({n[0]} {n[1]} {n[2]}) simpleGrading (1 1 1) );

boundary
(
{chr(10).join(pl)}
);
""")

    # --- system: the 1D suite's numerics, with GAMG and correctPhi for 3D AMR
    c = dict(maxCo=0.25, maxAlphaCo=0.25, maxAcousticCo=50, maxPhaseChangeFraction=0.05)
    c.update(spec.get("controls") or {})
    extra = "\n".join(f"{k:22s} {v};" for k, v in c.items())
    g.write(case, "system/controlDict", "dictionary", f"""application     compressibleLaserbeamFoam;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {spec['endTime']!r};
deltaT          1e-10;

writeControl    adjustableRunTime;
writeInterval   {spec['writeInterval']!r};
purgeWrite      0;
writeFormat     ascii;
writePrecision  12;
writeCompression off;
timeFormat      general;
timePrecision   10;
runTimeModifiable no;

adjustTimeStep  yes;
maxDeltaT       {spec['maxDeltaT']!r};
{extra}

// cell volumes and centres for the checker come from each written mesh (tools/check3d.py),
// serial or per processor, so AMR runs need no reconstructParMesh
""")
    scratch = os.path.join(OUT, ".scratch1d")
    g.OUT = scratch
    tmpl = g.make_case("tmpl", N=4, L=4e-6, gas=gas, left="wall", right="wall",
                       fields={"T": [300.0] * 4, "p": [1e5] * 4, "U": [0.0] * 4,
                               **g.phase_fields(4, gas)},
                       endTime=1.0, writeInterval=1.0, maxDeltaT=1.0,
                       liquidPr=spec.get("liquidPr", 0.16), vapourPr=spec.get("vapourPr", 0.67),
                       gasPr=spec.get("gasPr", 0.667))
    g.OUT = OUT
    with open(os.path.join(tmpl, "system/fvSchemes")) as f:
        schemes = f.read()
    with open(os.path.join(tmpl, "system/fvSolution")) as f:
        solution = f.read()
    solution = solution.replace("correctPhi               no;", "correctPhi               yes;")
    # GAMG in 3D. A near-uniform p or U (11, 12) has a round-off floor of ~1e-9 in OpenFOAM's
    # normalised residual, so the 1D tolerances would run every solve to maxIter
    for old, new in (("""    "pcorr.*"
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-10;
        relTol          0;
    }""", """    "pcorr.*"
    {
        solver          GAMG;
        smoother        DIC;
        tolerance       1e-9;
        relTol          0;
        maxIter         100;
    }"""), ("""    p_rgh
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-10;
        relTol          0;
    }""", """    p_rgh
    {
        solver          GAMG;
        smoother        DIC;
        tolerance       1e-9;
        relTol          0;
        maxIter         100;
    }"""), ("""        tolerance       1e-12;
        relTol          0;
        maxIter         200;""", """        tolerance       1e-11;
        relTol          0;
        maxIter         30;""")):
        assert solution.count(old) == 1, old
        solution = solution.replace(old, new)
    for rel, txt in (("system/fvSchemes", schemes), ("system/fvSolution", solution)):
        with open(os.path.join(case, rel), "w") as f:
            f.write(txt)
    g.write(case, "system/decomposeParDict", "dictionary", "numberOfSubdomains 2;\nmethod          scotch;\n")

    # --- constant: materials as in 1D; surface tension as asked
    sigma = spec.get("sigma", 0.0)
    os.makedirs(os.path.join(case, "constant"), exist_ok=True)
    for fn in os.listdir(os.path.join(tmpl, "constant")):
        src = os.path.join(tmpl, "constant", fn)
        if os.path.isfile(src) and fn.startswith("thermophysicalProperties"):
            with open(src) as f:
                txt = f.read()
            if fn == "thermophysicalProperties":
                txt = txt.replace("(metal1vapour metal1)   1.8", f"(metal1vapour metal1)   {sigma}")
                txt = txt.replace(f"({gas} metal1)          1.8", f"({gas} metal1)          {sigma}")
            with open(os.path.join(case, "constant", fn), "w") as f:
                f.write(txt)
    for fn in ("g", "turbulenceProperties", "LaserProperties", "timeVsLaserPosition", "timeVsLaserPower"):
        shutil.copy(os.path.join(tmpl, "constant", fn), os.path.join(case, "constant", fn))
    shutil.rmtree(scratch)

    a = dict(refineInterval=1, nBufferLayers=1, maxRefinement=2, maxCells=4000000, unrefineLevel=1e9)
    a.update(spec["amr"])
    g.write(case, "constant/dynamicMeshDict", "dictionary", f"""dynamicFvMesh   dynamicRefineFvMesh;

dynamicRefineFvMeshCoeffs
{{
    refineInterval      {a['refineInterval']};
    field               {a['field']};
    lowerRefineLevel    {a['lower']!r};
    upperRefineLevel    {a['upper']!r};
    unrefineLevel       {a['unrefineLevel']!r};
    nBufferLayers       {a['nBufferLayers']};
    maxRefinement       {a['maxRefinement']};
    maxCells            {a['maxCells']};
    // every flux the solver keeps over a mesh update; all are rebuilt from phi or the ledger
    correctFluxes
    (
        (phi none)
        (nHatf none)
        (rhoPhi none)
        (rhoCvPhi none)
        (rhoCvTPhi none)
        (ghf none)
        (alphaPhi0.metal1 none)
        (alphaPhi0.metal1vapour none)
        (alphaPhi0.{gas} none)
        (alphaPhi.metal1 none)
        (alphaPhi.metal1vapour none)
        (alphaPhi.{gas} none)
        (phiCorrDdt none)
        (ddtCoeff none)
        (phi_0 none)    // only written for restarts; on a dynamic mesh ddtCorr takes the old flux from Uf
        (uJump none)    // the Stefan jump per unit area (stefanFlow.H): mapped face by face, as stored
        (uJump_0 none)
        (phiJump none)  // recomputed every step before use
    );
    dumpLevel           false;
}}
""")

    # --- 0.orig/ on the base mesh (the Allrun copies it to 0/, then rewrites 0/ on the refined mesh)
    write_fields(case, spec, spec["rule"](centres(n, L), np.full(N, L[0] / n[0])), "0.orig")

    md = dict(name=name, dim=3, n=list(n), L=list(L), dx=L[0] / n[0], gas=gas,
              phases=["metal1", "metal1vapour", gas], maxRefinement=a["maxRefinement"], sigma=sigma,
              notes=spec.get("notes", ""))
    md.update(spec.get("meta") or {})
    with open(os.path.join(case, "case.json"), "w") as f:
        json.dump(md, f, indent=2)
    with open(os.path.join(case, "Allrun"), "w") as f:
        f.write(ALLRUN.replace("MAXREF", str(a["maxRefinement"])))
    os.chmod(os.path.join(case, "Allrun"), 0o755)
    return case


# ============================================================================= cases
def thermal_box3d():
    """hot argon sphere in cold argon, sealed octant, AMR on the thermal front"""
    Lb, R, Th, Tc = 100e-6, 50e-6, 1500.0, 300.0

    def rule(C, h):
        hot = fraction(C, h, lambda x, y, z: (x * x + y * y + z * z) < R * R)
        N = len(C)
        return {"T": 1.0 / (hot / Th + (1.0 - hot) / Tc),    # same p: the cell holds the mass of both parts
                "p": np.full(N, 1e5), "alpha.metal1": np.zeros(N), "alpha.metal1vapour": np.zeros(N),
                "alpha.argon": np.ones(N)}

    return dict(name="11_thermalBox3D", n=(16, 16, 16), L=(Lb, Lb, Lb), gas="argon", rule=rule,
                sides={"xmin": SYMMETRY, "ymin": SYMMETRY, "zmin": SYMMETRY,
                       "xmax": wall(), "ymax": wall(), "zmax": wall()},
                endTime=1.5e-5, writeInterval=5e-7, maxDeltaT=1e-7, gasPr=0.0667,
                amr=dict(field="T", lower=400.0, upper=1400.0),
                meta=dict(p0=1e5, Thot=Th, Tcold=Tc, R=R),
                notes="sealed adiabatic octant, hot argon sphere cooling by conduction while AMR refines "
                      "the thermal front and unrefines behind it: p must stay p0, mass and energy exact")


def contact3d():
    """hot Fe vapour sphere in cold argon (equal gamma) advected diagonally through a periodic box"""
    R, c0, U0, Th, Tc = 22e-6, (35e-6, 40e-6, 25e-6), (10.0, 5.0, 0.0), 3000.0, 300.0

    def rule(C, h):
        hot = fraction(C, h, lambda x, y, z: ((x - c0[0]) ** 2 + (y - c0[1]) ** 2 + (z - c0[2]) ** 2) < R * R)
        N = len(C)
        # each part keeps its own volume at the common p, and the mixed T keeps that p
        return {"T": 1.0 / (hot / Th + (1.0 - hot) / Tc), "p": np.full(N, 1e5), "U": U0,
                "alpha.metal1": np.zeros(N), "alpha.metal1vapour": hot, "alpha.argon": 1.0 - hot}

    return dict(name="12_contact3D", n=(16, 16, 8), L=(100e-6, 100e-6, 50e-6), gas="argon", rule=rule,
                sides={"xmin": CYCLIC, "xmax": CYCLIC, "ymin": CYCLIC, "ymax": CYCLIC,
                       "zmin": SYMMETRY, "zmax": SYMMETRY},
                endTime=1e-5, writeInterval=1e-6, maxDeltaT=1e-7, gasPr=1e6, vapourPr=1e6,
                amr=dict(field="T", lower=350.0, upper=2950.0),
                meta=dict(p0=1e5, U0=list(U0), Thot=Th, Tcold=Tc, R=R, c0=list(c0)),
                notes="hot Fe vapour sphere in cold argon (both gamma 5/3) carried at U0 = (10, 5, 0) m/s "
                      "through a box periodic in x and y, conduction off, AMR following the sphere: U and p "
                      "stay uniform, momentum and mass exact, the sphere moves at U0")


def evap_closed3d(kind):
    """sealed box, liquid Fe under pure Fe vapour away from saturation, walls at T0"""
    T0, p0, Lb = 3100.0, (0.93e5 if kind == "dropletCond" else 0.8e5), 100e-6
    if kind == "film":
        z0, a, b = 40e-6, 0.4, 0.2
        inside = lambda x, y, z: z < z0 + a * (x - 0.5 * Lb) + b * (y - 0.5 * Lb)   # noqa: E731
        sides = {s: wall(T0) for s in SIDES}
        area = Lb * Lb * math.sqrt(1 + a * a + b * b)
        name, geo, endTime, wint = "13_inclinedFilm3D", dict(z0=z0, a=a, b=b), 3e-6, 2e-7
        notes = ("liquid Fe film whose surface is the plane z = z0 + a(x - L/2) + b(y - L/2), inclined to the "
                 "mesh, under pure Fe vapour at 0.8 bar in a sealed box at 3100 K: the evaporation rate must "
                 "be j*A for the true plane area, the end state Psat(3100 K)")
    else:
        R = 40e-6
        inside = lambda x, y, z: (x * x + y * y + z * z) < R * R                   # noqa: E731
        sides = {"xmin": SYMMETRY, "ymin": SYMMETRY, "zmin": SYMMETRY,
                 "xmax": wall(T0), "ymax": wall(T0), "zmax": wall(T0)}
        area = 4 * math.pi * R * R / 8
        geo, endTime, wint = dict(R=R), 1e-5, 5e-7
        if kind == "dropletCond":
            name = "14_dropletCond3D"
            notes = ("octant of a liquid Fe droplet (R = 40 um) in supersaturated pure Fe vapour at 0.93 bar, sealed "
                     "box with walls at 3100 K: condensation on every interface orientation and AMR level at "
                     "j*4 pi R^2/8, end state Psat(T_liquid), droplet stays round and still")
        else:
            name = "16_dropletEvap3D"
            notes = ("14 run the other way, evaporating at 0.8 bar: the droplet is superheated (p < Psat(T)) until the "
                     "box saturates, and with no nucleation barrier vapour traces inside it grow (known limitation)")

    def rule(C, h):
        cond = fraction(C, h, inside)
        N = len(C)
        return {"T": np.full(N, T0), "p": np.full(N, p0), "alpha.metal1": cond, "alpha.metal1vapour": 1.0 - cond,
                "alpha.argon": np.zeros(N)}

    return dict(name=name, n=(16, 16, 16), L=(Lb, Lb, Lb), gas="argon", rule=rule, sides=sides,
                endTime=endTime, writeInterval=wint, maxDeltaT=1e-7, liquidPr=0.016, vapourPr=0.0067,
                amr=dict(field="alpha.metal1", lower=0.001, upper=0.999, unrefineLevel=10),
                meta=dict(T0=T0, p0=p0, Psat0=g.psat(T0), area=area, kind=kind, **geo),
                notes=notes + " (surface tension off so the interface stays where it is; conductivities "
                              "raised as in 07)")


def evap_steady3d():
    """test 09 in 3D: heated film, open vapour column, AMR on the interface"""
    h_l, q = 100e-6, 1e8
    kl = 6.0e-3 * 800 / 0.016
    j = q / 6.09e6
    lo, hi = 3000.0, 5000.0
    for _ in range(200):
        Ts0 = 0.5 * (lo + hi)
        rv = 1e5 / (8314.47 / g.M_FE * Ts0)
        jh = math.sqrt((g.M_FE / 1000) / (2 * math.pi * g.RU * Ts0)) * (g.psat(Ts0) - 1e5 - j * j / rv)
        lo, hi = (Ts0, hi) if jh < j else (lo, Ts0)
    Ts0 = round(Ts0, 3)

    def rule(C, h):
        z = C[:, 2]
        cond = np.clip((h_l - (z - 0.5 * h)) / h, 0.0, 1.0)
        N = len(z)
        return {"T": Ts0 + q * np.maximum(h_l - z, 0.0) / kl,    # steady conduction profile in the film
                "p": np.full(N, 1e5), "alpha.metal1": cond, "alpha.metal1vapour": 1.0 - cond,
                "alpha.argon": np.zeros(N)}

    bottom = {"_type": "wall",
              "U": "type fixedValue; value uniform (0 0 0);",
              "p_rgh": "type fixedFluxPressure; value $internalField;",
              "p": "type calculated; value $internalField;",
              "T": f"type fixedGradient; gradient uniform {q / kl!r};",
              "Laser_boundary": "type fixedValue; value uniform -1;"}
    top = {"_type": "patch",
           "U": "type pressureInletOutletVelocity; value uniform (0 0 0);",
           "p_rgh": "type fixedValue; value uniform 1e5;",
           "p": "type calculated; value $internalField;",
           "T": f"type inletOutlet; inletValue uniform {Ts0}; value uniform {Ts0};",
           "alpha.metal1": "type inletOutlet; inletValue uniform 0; value uniform 0;",
           "alpha.metal1vapour": "type inletOutlet; inletValue uniform 1; value uniform 1;",
           "alpha.argon": "type inletOutlet; inletValue uniform 0; value uniform 0;",
           "Laser_boundary": "type fixedValue; value uniform 1;"}
    return dict(name="15_evapSteady3D", n=(4, 4, 100), L=(40e-6, 40e-6, 1e-3), gas="argon", rule=rule,
                sides={"xmin": SYMMETRY, "xmax": SYMMETRY, "ymin": SYMMETRY, "ymax": SYMMETRY,
                       "zmin": bottom, "zmax": top},
                endTime=3e-5, writeInterval=2.5e-6, maxDeltaT=1e-7, liquidPr=0.016, sigma=1.8,
                amr=dict(field="alpha.metal1", lower=0.001, upper=0.999, unrefineLevel=10, maxRefinement=1),
                meta=dict(q=q, kl=kl, h_l=h_l, p_inf=1e5, Ts0=Ts0),
                notes="test 09 on a 3D column: liquid Fe film heated from below at 1e8 W/m2, vapour column "
                      "open at 1 bar, AMR on the interface (5 um there, 10 um in the column): jump conditions. "
                      "Surface tension on (1.8 N/m, as in 09): the film is 60-90 K superheated, and with sigma = 0 "
                      "it boils from grid-scale perturbations within 1 us in 3D")


SPECS = (thermal_box3d, contact3d, lambda: evap_closed3d("film"), lambda: evap_closed3d("dropletCond"), evap_steady3d,
         lambda: evap_closed3d("dropletEvap"))


# ============================================================================= --orig, --init
def spec_of(case):
    try:
        name = json.load(open(os.path.join(case, "case.json")))["name"]
    except (OSError, ValueError, KeyError):
        name = os.path.basename(os.path.normpath(case))      # case.json lost too: go by the folder name
    return next(s for s in (f() for f in SPECS) if s["name"] == name)


def orig(case):
    """Rewrite 0.orig, the fields on the base mesh, from the case's spec."""
    spec = spec_of(case)
    n, L = spec["n"], spec["L"]
    write_fields(case, spec, spec["rule"](centres(n, L), np.full(n[0] * n[1] * n[2], L[0] / n[0])), "0.orig")
    print(f"{spec['name']}: wrote 0.orig on the {n[0]}x{n[1]}x{n[2]} base mesh")


def init(case):
    """Move the mesh the solver refined to t = 0 and write the fields on it."""
    import foamio
    import mesh3d

    spec = spec_of(case)
    name = spec["name"]
    ts = [t for t in foamio.times(case) if float(t) > 0]
    if not ts:
        sys.exit("the refining steps wrote nothing (see log.prerefine)")
    src, dst = os.path.join(case, ts[-1], "polyMesh"), os.path.join(case, "constant", "polyMesh")
    if os.path.isfile(os.path.join(src, "faces")):              # no polyMesh there: nothing was refined
        for fn in os.listdir(src):
            shutil.copy(os.path.join(src, fn), os.path.join(dst, fn))
    for t in ts:
        shutil.rmtree(os.path.join(case, t))
    V, C = mesh3d.mesh(case, "0")
    h = np.cbrt(V)
    shutil.rmtree(os.path.join(case, "0"))          # also the fields the refining steps wrote at t = 0
    write_fields(case, spec, spec["rule"](C, h))
    lv = np.round(np.log2(spec["L"][0] / spec["n"][0] / h)).astype(int)
    print(f"{name}: {len(V)} cells at t = 0, per level {np.bincount(lv).tolist()}, "
          f"volume {V.sum():.6e} (block {np.prod(spec['L']):.6e})")


if __name__ == "__main__":
    if MODE == "--init":
        init(os.path.abspath(TARGET))
        sys.exit(0)
    if MODE == "--orig":
        orig(os.path.abspath(TARGET))
        sys.exit(0)
    if MODE == "--setup":
        spec = spec_of(os.path.abspath(TARGET))
        OUT = os.path.dirname(os.path.abspath(os.path.normpath(TARGET)))
        print("wrote", make_case3d(spec))
        sys.exit(0)
    os.makedirs(OUT, exist_ok=True)
    for s in SPECS:
        print("wrote", make_case3d(s()))
