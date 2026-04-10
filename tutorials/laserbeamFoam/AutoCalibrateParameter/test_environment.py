#!/usr/bin/env python3
"""
Environment verification for the AutoCalibrateParameter tutorial.

This script checks that all required Python packages are installed, the
calibration framework imports correctly, and — crucially — runs the meltpool
post-processing pipeline end-to-end on bundled test data to verify that the
depth/width extraction chain (numpy → pandas → joblib → matplotlib) works.

Usage:
    python test_environment.py
"""
import importlib
import os
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

# ---------------------------------------------------------------------------
# Resolve paths
# ---------------------------------------------------------------------------
TUTORIAL_DIR = Path(__file__).resolve().parent
REPO_ROOT = TUTORIAL_DIR.parents[2]                          # -> LaserbeamFoam
SCRIPTS_DIR = REPO_ROOT / "applications" / "scripts"
PKG_ROOT = SCRIPTS_DIR / "AutoCalibrateParameter"
POSTPROC_DIR = REPO_ROOT / "applications" / "scripts" / "postProcessing"

# Make AutoCalibrateParameter importable
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

# ---------------------------------------------------------------------------
# Counters
# ---------------------------------------------------------------------------
_pass = 0
_fail = 0
_warn = 0


def ok(msg):
    global _pass
    _pass += 1
    print(f"  [PASS] {msg}")


def fail(msg):
    global _fail
    _fail += 1
    print(f"  [FAIL] {msg}")


def warn(msg):
    global _warn
    _warn += 1
    print(f"  [WARN] {msg}")


# ===================================================================
# 1. Required Python packages
# ===================================================================
def test_python_packages():
    print("\n1. Checking required Python packages …")
    required = {
        "numpy": "array computation",
        "pandas": "CSV / data-frame handling",
        "yaml": "YAML config parsing (PyYAML)",
        "matplotlib": "plotting convergence curves",
        "joblib": "meltpool geometry caching (used by postprocessing)",
        "scipy": "gradient-based optimisation",
    }
    optional = {
        "skopt": "Bayesian optimisation (scikit-optimize)",
    }

    for mod, purpose in required.items():
        try:
            importlib.import_module(mod)
            ok(f"{mod} — {purpose}")
        except ImportError:
            fail(f"{mod} — {purpose}  ← pip install {mod}")

    for mod, purpose in optional.items():
        try:
            importlib.import_module(mod)
            ok(f"{mod} — {purpose}")
        except ImportError:
            warn(f"{mod} — {purpose} (needed only for --method bayes)")


# ===================================================================
# 2. AutoCalibrateParameter package import
# ===================================================================
def test_package_import():
    print("\n2. Importing AutoCalibrateParameter …")
    try:
        from AutoCalibrateParameter import (
            BaseConfig,
            BayesConfig,
            GradientConfig,
            BayesOptimizer,
            GradientOptimizer,
            SimulationRunner,
            load_experiment_data,
        )
        ok("All core classes imported successfully")
    except ImportError as e:
        fail(f"Import error: {e}")


# ===================================================================
# 3. Config parsing
# ===================================================================
def test_config_parsing():
    print("\n3. Parsing config.yaml …")
    config_path = TUTORIAL_DIR / "config.yaml"
    if not config_path.exists():
        fail(f"config.yaml not found at {config_path}")
        return

    try:
        from AutoCalibrateParameter.config import BayesConfig
        config = BayesConfig.from_yaml(config_path)
        ok(f"Config loaded: {len(config.active_params)} active params, "
           f"n_batches={config.n_batches}")
    except Exception as e:
        fail(f"Config parsing failed: {e}")


# ===================================================================
# 4. Experimental data loading
# ===================================================================
def test_experimental_data():
    print("\n4. Loading experimental data …")
    csv_path = TUTORIAL_DIR / "SingleTrackExperimentalData.csv"
    if not csv_path.exists():
        fail(f"CSV not found: {csv_path}")
        return

    try:
        import pandas as pd
        df = pd.read_csv(csv_path)
        required_cols = {"power_W", "depth_um", "width_um"}
        missing = required_cols - set(df.columns)
        if missing:
            fail(f"Missing columns: {missing}")
        else:
            ok(f"Loaded {len(df)} power conditions: "
               f"{', '.join(str(int(p)) + 'W' for p in df['power_W'])}")
    except Exception as e:
        fail(f"CSV load error: {e}")


# ===================================================================
# 5. Post-processing dependency chain
# ===================================================================
def test_postproc_dependencies():
    print("\n5. Testing meltpool post-processing dependencies …")

    postproc_script = POSTPROC_DIR / "characterise_meltpool.py"
    if postproc_script.exists():
        ok("Post-processing script found")
    else:
        fail(f"Post-processing script missing: {postproc_script}")

    functions_py = POSTPROC_DIR / "functions.py"
    if functions_py.exists():
        ok("functions.py found")
    else:
        fail(f"functions.py missing: {functions_py}")

    # Test the critical chain: numpy + pandas + joblib + matplotlib
    try:
        import numpy as np
        import pandas as pd
        from joblib import dump, load
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("x,y,z,alpha.metal\n")
            f.write("0.00015,0.0004,0.00015,1.0\n")
            f.write("0.00015,0.0004,0.00014,0.8\n")
            f.write("0.00015,0.0004,0.00013,0.0\n")
            tmp_csv = f.name

        df = pd.read_csv(tmp_csv)
        arr = df[["x", "y", "z", "alpha.metal"]].to_numpy()
        assert arr.shape == (3, 4)

        with tempfile.NamedTemporaryFile(suffix=".joblib", delete=False) as jf:
            dump(arr, jf.name)
            loaded = load(jf.name)
            assert np.allclose(arr, loaded)
            os.unlink(jf.name)

        os.unlink(tmp_csv)
        ok("numpy + pandas + joblib + matplotlib chain works")
    except Exception as e:
        fail(f"Post-processing dependency chain broken: {e}")

    if shutil.which("pvpython"):
        ok("pvpython found in PATH")
    else:
        warn("pvpython not in PATH — needed for VTK-to-CSV extraction. "
             "Install ParaView or add pvpython to PATH.")


# ===================================================================
# 6. OpenFOAM case integrity
# ===================================================================
def test_openfoam_case():
    print("\n6. Checking OpenFOAM case files …")
    required = {
        "initial/alpha.metal": "initial volume fraction field",
        "initial/T": "initial temperature field",
        "initial/U": "initial velocity field",
        "initial/p_rgh": "initial pressure field",
        "constant/transportProperties": "material properties",
        "constant/LaserProperties": "laser configuration",
        "constant/location": "powder bed particle positions",
        "constant/timeVsLaserPower": "laser power schedule",
        "constant/timeVsLaserPosition": "laser scan path",
        "system/controlDict": "solver control",
        "system/blockMeshDict": "mesh definition",
        "system/fvSchemes": "discretisation schemes",
        "system/fvSolution": "solver settings",
        "system/decomposeParDict": "parallel decomposition",
    }
    for rel, desc in required.items():
        p = TUTORIAL_DIR / rel
        if p.exists():
            ok(f"{rel}")
        else:
            fail(f"{rel} — {desc}")


# ===================================================================
# 7. End-to-end post-processing on real simulation results
# ===================================================================
def test_postproc_end_to_end():
    """Run characterise_meltpool.py on bundled simulation output and compare.

    test_data/ contains:
      - 0.0007/alpha.metal  — a real OpenFOAM simulation result
      - constant/polyMesh/  — the mesh definition
      - system/controlDict  — minimal solver dict (for OpenFOAMReader)
      - main.foam           — empty marker file for ParaView
      - expected_metrics_summary.csv — baseline values

    The test runs the FULL pipeline:
      pvpython (VTK extraction) → geometry analysis → baseline comparison

    If pvpython is not available the test is skipped with a warning.
    """
    print("\n7. Running post-processing on test data (end-to-end) …")

    test_data_dir = TUTORIAL_DIR / "test_data"
    if not test_data_dir.exists():
        warn("test_data/ directory not found — skipping end-to-end test")
        return

    expected_csv = test_data_dir / "expected_metrics_summary.csv"
    if not expected_csv.exists():
        fail("test_data/expected_metrics_summary.csv missing")
        return

    # Check that simulation result exists
    alpha_field = test_data_dir / "0.0007" / "alpha.metal"
    if not alpha_field.exists():
        fail("test_data/0.0007/alpha.metal missing — no simulation result to test")
        return

    if not shutil.which("pvpython"):
        warn("pvpython not available — skipping end-to-end extraction test")
        return

    # Generate input_data.py required by characterise_meltpool.py
    input_data_path = test_data_dir / "input_data.py"
    input_data_path.write_text(textwrap.dedent("""\
        LASER_DIAMETER = 7e-05
        OF_LOCATION = ''
        CELL_SIZE = 1e-06
        X_MIN_AND_MAX_DOMAIN = [0.0, 0.0003]
        Y_COORD_BEGIN_TRACK = 0.00015
        Y_COORD_END_TRACK = 0.0003
        PLOT_GEOMETRY_VS_Y_LOCATION = False
    """))

    postproc_script = POSTPROC_DIR / "characterise_meltpool.py"

    # Clean previous outputs so we get a fresh full run
    for name in ("meltpool.csv", "meltpool_slice_xmid.csv",
                 "cross_sections_statistics.csv", "row_statistics.csv",
                 "metrics_summary.csv", "metrics_summary.xlsx",
                 "continuous.joblib"):
        p = test_data_dir / name
        if p.exists():
            p.unlink()

    env = os.environ.copy()
    # Do NOT set SKIP_PVPYTHON — run the full pipeline including VTK extraction

    try:
        result = subprocess.run(
            [sys.executable, str(postproc_script)],
            cwd=str(test_data_dir),
            env=env,
            capture_output=True,
            text=True,
            timeout=120,
        )
        if result.returncode != 0:
            fail(f"characterise_meltpool.py exited with code {result.returncode}")
            if result.stderr:
                print(f"        stderr: {result.stderr.strip()[:500]}")
            return
    except subprocess.TimeoutExpired:
        fail("Post-processing timed out (>120 s)")
        return
    except Exception as e:
        fail(f"Failed to run post-processing: {e}")
        return
    finally:
        # Clean up generated files
        if input_data_path.exists():
            input_data_path.unlink()

    # Verify meltpool.csv was extracted by pvpython
    meltpool_csv = test_data_dir / "meltpool.csv"
    if not meltpool_csv.exists():
        fail("pvpython did not generate meltpool.csv from simulation data")
        return
    ok("pvpython extracted meltpool.csv from simulation result")

    # Verify metrics output
    output_csv = test_data_dir / "metrics_summary.csv"
    if not output_csv.exists():
        fail("metrics_summary.csv was not generated")
        return
    ok("Geometry analysis completed → metrics_summary.csv")

    # Compare against baseline
    try:
        import pandas as pd
        df_expected = pd.read_csv(expected_csv)
        df_actual = pd.read_csv(output_csv)

        cols = ["width_mean_um", "height_mean_um", "depth_mean_um"]
        all_match = True
        for col in cols:
            val_exp = df_expected[col].iloc[0]
            val_act = df_actual[col].iloc[0]
            diff = abs(val_exp - val_act)
            if diff > 0.5:  # 0.5 µm tolerance
                fail(f"{col}: expected {val_exp:.2f}, got {val_act:.2f} "
                     f"(Δ {diff:.2f} µm)")
                all_match = False
            else:
                ok(f"{col}: {val_act:.2f} µm ✓")

        if all_match:
            ok("All metrics match baseline — full pipeline works!")
    except Exception as e:
        fail(f"Baseline comparison failed: {e}")


# ===================================================================
# Summary
# ===================================================================
def main():
    print("=" * 64)
    print("  AutoCalibrateParameter — Environment Verification")
    print("=" * 64)
    print(f"  Python:    {sys.version.split()[0]}")
    print(f"  Repo root: {REPO_ROOT}")
    print(f"  Tutorial:  {TUTORIAL_DIR}")

    test_python_packages()
    test_package_import()
    test_config_parsing()
    test_experimental_data()
    test_postproc_dependencies()
    test_openfoam_case()
    test_postproc_end_to_end()

    print()
    print("=" * 64)
    if _fail == 0:
        tag = "ALL CHECKS PASSED" + (f" ({_warn} warnings)" if _warn else "")
        print(f"  ✓ {tag}")
        print("  You are ready to run:  python main.py")
    else:
        print(f"  ✗ {_fail} FAILED, {_pass} passed, {_warn} warnings")
        print("  Please fix the failures above before running the calibration.")
    print("=" * 64)

    return 1 if _fail else 0


if __name__ == "__main__":
    sys.exit(main())
