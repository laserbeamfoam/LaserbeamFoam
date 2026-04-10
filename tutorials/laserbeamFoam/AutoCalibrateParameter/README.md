# AutoCalibrateParameter — Meltpool Calibration Tutorial

This tutorial demonstrates how to automatically calibrate `laserbeamFoam`
simulation parameters against experimental meltpool measurements using the
**AutoCalibrateParameter** framework.

The case models a **316L stainless steel single laser track** at 200 W.
The optimiser iteratively runs OpenFOAM simulations, extracts the meltpool
width and depth via the post-processing pipeline, and adjusts model parameters
to minimise the discrepancy with experimental data.

---

## How It Works

In laser powder bed fusion (LPBF) simulations, several physical parameters
(surface tension, Marangoni coefficient, recoil pressure coefficient, etc.)
are difficult to determine from first principles. This framework automates
their calibration by:

1. **Propose** — the optimiser (Bayesian or gradient-based, see below)
   suggests a set of candidate parameter values (e.g. absorptivity, Marangoni
   coefficient, surface tension, recoil pressure coefficient).
2. **Simulate** — run `laserbeamFoam` via `mpirun` for each experimental
   power setting (200 W). Mesh generation and field initialisation are handled
   automatically by the framework.
3. **Post-process** — call `characterise_meltpool.py` (via `pvpython`) to
   extract meltpool geometry from the simulation results.
4. **Evaluate** — compute the NRMSE between simulated and experimental
   width/depth.
5. **Update** — feed the result back into the optimiser and repeat until
   convergence.

### Optimisation algorithms

| Method | Flag | Description |
|---|---|---|
| **Bayesian** (default) | `--method bayes` | Gaussian process surrogate + acquisition function (EI/PI/LCB). Efficient for expensive black-box functions; builds a probabilistic model of the objective and balances exploration vs. exploitation. |
| **Gradient-based** | `--method gradient` | `scipy.optimize.least_squares` with finite-difference gradients. Faster per iteration but may converge to local minima; supports multi-start to mitigate this. |

Calibrated parameter values, convergence history, and diagnostic plots are
saved under `runs/`.

---

## Directory Structure

```
AutoCalibrateParameter/
├── Allinstall                         # One-shot install script — run this first
├── main.py                            # Entry point — run this
├── config.yaml                        # Calibration settings (editable)
├── test_environment.py                # Environment check
├── SingleTrackExperimentalData.csv    # Target experimental data
├── Allrun                             # OpenFOAM case setup script
├── Allclean                           # OpenFOAM case cleanup script
├── initial/                           # Initial boundary conditions
├── constant/                          # Material properties & laser config
│   ├── transportProperties            #   SS316L thermo-physical properties
│   ├── LaserProperties                #   Laser model parameters
│   ├── location                       #   Powder bed particle positions
│   ├── timeVsLaserPower               #   Laser power schedule
│   └── timeVsLaserPosition            #   Laser scan path
└── system/                            # Solver settings
    ├── controlDict                    #   Time stepping & output control
    ├── blockMeshDict                  #   Mesh definition (300×800×300 µm)
    ├── decomposeParDict               #   Parallel decomposition
    ├── fvSchemes
    ├── fvSolution
    ├── setFieldsDict
    └── bedPlateDict
```

---

## Full Installation Guide

Follow these steps in order before running the calibration for the first time.

### Step 1 — Install ParaView (with pvpython)

`pvpython` is required for meltpool post-processing. The recommended way is
via conda, which also bundles `pvpython` automatically:

```bash
conda install -c conda-forge paraview=6.0.1
```

After installation, verify:

```bash
pvpython --version
```

> If you install ParaView from the system package manager or official
> installer instead, make sure `pvpython` is on your `PATH`.

### Step 2 — Create the Python environment

Create a dedicated conda environment with all required packages pinned to
the tested versions:

```bash
conda create -n meltpool-postproc python=3.10.19 \
    -c conda-forge \
    numpy=2.2.6 \
    pandas=2.3.3 \
    matplotlib=3.10.8 \
    joblib=1.5.2 \
    paraview=6.0.1

conda activate meltpool-postproc

pip install \
    scipy==1.15.3 \
    PyYAML==6.0.3 \
    openpyxl==3.1.5 \
    scikit-optimize==0.10.2
```

> `scikit-optimize` is only required for `--method bayes` (the default).
> For `--method gradient`, you can omit it.

### Step 3 — Source your OpenFOAM environment

```bash
of2506    # or: source $WM_PROJECT_DIR/etc/bashrc
```

### Step 4 — Run Allinstall

With the conda environment active and OpenFOAM sourced, run:

```bash
conda activate meltpool-postproc
./Allinstall
```

This will:
1. Install the `AutoCalibrateParameter` package into the active environment
   (editable install — source changes take effect immediately)
2. Copy the `postProcessing` scripts to `$FOAM_USER_APPBIN`
3. Verify all imports pass

### Step 5 — Verify

```bash
python test_environment.py
```

This checks Python packages, external tools (`pvpython`, `blockMesh`,
`laserbeamFoam`, `mpirun`), and the `AutoCalibrateParameter` package imports.

---

## Quick Start (after installation)

```bash
# Activate environment and source OpenFOAM
conda activate meltpool-postproc
of2506

# Run calibration (Bayesian, default)
python main.py

# Or gradient-based
python main.py --method gradient

# See all options
python main.py --help
```

> **No manual case setup needed.** `main.py` calls `blockMesh`,
> `setSolidFraction`, and `cp -r initial 0` automatically before each
> simulation run.

### Production runs

The default `config.yaml` is tuned for a **quick demonstration** (5 initial
points, 5 batches). For production-quality calibration, edit `config.yaml`:

```yaml
n_initial_points: 20   # more initial exploration
n_batches: 30          # more optimisation iterations
n_proc: 24             # match your CPU count
```

---

## Common Options

| Flag | Default | Description |
|---|---|---|
| `--method {bayes,gradient}` | `bayes` | Optimisation algorithm |
| `--config FILE` | `config.yaml` | Path to configuration file |
| `--n-proc N` | 4 | Number of MPI ranks |
| `--n-initial-points N` | 5 | Random initial samples (Bayes) |
| `--n-batches N` | 5 | Optimisation iterations (Bayes) |
| `--hpc` | off | Enable HPC job submission mode |
| `--verbose` | off | Detailed logging |

---

## Configuration Reference

All calibration parameters are set in `config.yaml`:

| Key | Default | Description |
|---|---|---|
| `n_initial_points` | 5 | Initial random samples (Bayesian) |
| `n_batches` | 5 | Optimisation iterations |
| `n_proc` | 4 | MPI processes for OpenFOAM |
| `active_params` | sigma, marangoni, … | Parameters to optimise |
| `fixed_values` | substrate_temp: 300 | Held-constant parameters |
| `output_weights` | [1, 1, 0] | Weights for [width, depth, area] in cost function |
| `postproc_script` | (auto-resolved) | Path to meltpool extraction script |

- **`active_params`** — choose which parameters to optimise. Inactive
  parameters are held at the values in `fixed_values` (or midpoint of their
  bounds).
- **`output_weights`** — set a component to `0` to exclude it from the cost
  function (e.g. `[1, 1, 0]` ignores meltpool area, which is the default).

---

## Experimental Data Format

`SingleTrackExperimentalData.csv` must contain at least these columns:

```csv
power_W,depth_um,width_um
200,56.87,111.52
```

Each row represents one laser power condition with its measured meltpool
depth and width in micrometres. The optimiser runs a simulation for each row
and compares the predicted geometry against the measured values. The tutorial
uses a single condition at 200 W; add more rows to calibrate across multiple
power levels.

---

## Troubleshooting

| Symptom | Solution |
|---|---|
| `pvpython: command not found` | Install ParaView via conda (`conda install -c conda-forge paraview`) and ensure the env is active |
| `blockMesh not found` | Source the OpenFOAM environment (`of2506` or equivalent) |
| `ImportError: No module named 'skopt'` | `pip install scikit-optimize==0.10.2` |
| `ImportError: No module named 'yaml'` | `pip install PyYAML==6.0.3` |
| `ImportError: No module named 'openpyxl'` | `pip install openpyxl==3.1.5` |
| `meltpool.csv` is empty | Ensure simulation ran to completion |
| Post-processing returns NaN | Check `y_begin_track` / `y_end_track` in config are inside the domain |
| `foam_runner` error | Set `foam_runner: null` if OpenFOAM is already sourced |
| Post-processing fails silently | Run `python test_environment.py` for diagnostics |

---

## References

- [LaserbeamFoam repository](https://github.com/laserbeamfoam/LaserbeamFoam)
- [scikit-optimize documentation](https://scikit-optimize.github.io/)
- Flint, T.F. et al., *A fundamental investigation into the role of beam
  focal point, and beam divergence, on thermo-capillary stability and
  evolution of electron beam welding induced porosity*, International Journal
  of Heat and Mass Transfer, 2021.
