# AutoCalibrateParameter — Meltpool Calibration Tutorial

This tutorial demonstrates how to automatically calibrate `laserbeamFoam`
simulation parameters against experimental meltpool measurements using the
**AutoCalibrateParameter** framework.

The case is a **316L stainless steel single track** laser scan setup. The
optimiser iteratively runs OpenFOAM simulations, extracts the meltpool width
and depth via the post-processing pipeline, and adjusts model parameters to
minimise the discrepancy with experimental data.

---

## Prerequisites

| Requirement | Version | Notes |
|---|---|---|
| OpenFOAM | v2412 / v2506 | `laserbeamFoam` solver must be compiled |
| ParaView | ≥ 5.7 | `pvpython` must be in `PATH` |
| Python | ≥ 3.10 | with packages listed below |
| MPI | any | `mpirun` for parallel runs |

### Python packages

```bash
# Option A: use the provided conda environment
conda env create -f ../../../applications/scripts/postProcessing/environment.yml
conda activate meltpool-postproc
pip install scikit-optimize scipy pyyaml

# Option B: install into an existing environment
pip install numpy pandas matplotlib joblib scikit-optimize scipy pyyaml
```

### Verify the environment

Before running the calibration, check that everything is in place:

```bash
python test_environment.py
```

This script verifies Python packages, external tools (`pvpython`,
`blockMesh`, `laserbeamFoam`, `mpirun`), the post-processing script, and the
`AutoCalibrateParameter` package imports.

---

## Directory Structure

```
AutoCalibrateParameter/
├── Allrun                             # OpenFOAM case setup script
├── Allclean                           # OpenFOAM case cleanup script
├── config.yaml                        # Calibration configuration
├── main.py                            # Entry point
├── test_environment.py                # Environment check
├── input_data.py                      # Post-processing geometry settings
├── SingleTrackExperimentalData.csv    # Experimental measurements
├── README.md                          # This file
├── initial/                           # Initial field files
├── constant/                          # Physical properties & laser settings
│   ├── LaserProperties
│   ├── transportProperties
│   ├── timeVsLaserPosition
│   └── timeVsLaserPower
└── system/                            # Solver & mesh settings
    ├── blockMeshDict
    ├── controlDict
    ├── fvSchemes
    ├── fvSolution
    ├── decomposeParDict
    ├── setFieldsDict
    └── bedPlateDict
```

---

## How It Works

1. **Experimental data** — Three single-track experiments at different laser
   powers (140 W, 200 W, 260 W) provide target meltpool width and depth
   values.

2. **Optimiser loop** — At each iteration the framework:
   - Proposes a set of simulation parameters (e.g. absorption coefficient,
     Marangoni coefficient, surface tension, recoil pressure coefficient).
   - Runs the `laserbeamFoam` solver via `mpirun` for each experimental
     power setting.
   - Calls `characterise_meltpool.py` (with `pvpython`) to extract meltpool
     geometry from the simulation results.
   - Computes the NRMSE between simulated and experimental width/depth.
   - Updates the parameter estimates.

3. **Output** — Calibrated parameter values, convergence history, and
   diagnostic plots are saved under `runs/`.

---

## Quick Start

```bash
# 1. Set up the case
source $WM_PROJECT_DIR/bin/tools/RunFunctions
./Allrun

# 2. Run the calibration (Bayesian optimisation, default)
python main.py

# 3. Or use gradient-based optimisation
python main.py --method gradient
```

### Common options

| Flag | Default | Description |
|---|---|---|
| `--method {bayes,gradient}` | `bayes` | Optimisation algorithm |
| `--config FILE` | `config.yaml` | Path to configuration file |
| `--n-proc N` | 24 | Number of MPI ranks |
| `--n-initial-points N` | 20 | Random initial samples (Bayes) |
| `--n-batches N` | 20 | Optimisation iterations (Bayes) |
| `--hpc` | off | Enable HPC job submission mode |
| `--verbose` | off | Detailed logging |

See all options with `python main.py --help`.

---

## Configuration

Edit `config.yaml` to customise the calibration. Key settings:

- **`active_params`** — choose which parameters to optimise. Inactive
  parameters are held at the values in `fixed_values` (or midpoint of their
  bounds).
- **`output_weights`** — weight for `[width, depth, area]` in the cost
  function. Set a component to `0` to exclude it.
- **`n_proc`** — adjust to match your machine's core count.
- **`postproc_script`** — relative path to `characterise_meltpool.py`.

---

## Experimental Data Format

The CSV file must contain at least these columns:

```csv
power_W,depth_um,width_um
140,47.53,89.61
200,56.87,111.52
260,63.15,133.1
```

Each row is one experimental condition. The optimiser runs a separate
simulation for each power setting and compares the predicted meltpool
geometry against the measured values.

---

## Troubleshooting

| Problem | Solution |
|---|---|
| `pvpython not found` | Install ParaView or add it to `PATH` |
| `blockMesh not found` | Source OpenFOAM environment (`of2506` or equivalent) |
| `ImportError: No module named 'skopt'` | `pip install scikit-optimize` |
| `ImportError: No module named 'yaml'` | `pip install pyyaml` |
| Post-processing fails silently | Run `python test_environment.py` for diagnostics |

---

## References

- [laserbeamFoam on GitHub](https://github.com/laserbeamfoam/LaserbeamFoam)
- Flint, T.F. et al., *A fundamental investigation into the role of beam
  focal point, and beam divergence, on thermo-capillary stability and
  evolution of electron beam welding induced porosity*, International Journal
  of Heat and Mass Transfer, 2021.
# AutoCalibrateParameter Tutorial — 316L Single Track Calibration

This tutorial demonstrates how to automatically calibrate `laserbeamFoam`
simulation parameters against experimental meltpool measurements using the
**AutoCalibrateParameter** framework.

## What This Does

In laser powder bed fusion (LPBF) simulations, several physical parameters
(surface tension, Marangoni coefficient, recoil pressure coefficient, etc.)
are difficult to determine from first principles. This framework automates
their calibration by:

1. Running `laserbeamFoam` with candidate parameter sets
2. Post-processing each simulation to extract predicted meltpool width/depth
3. Comparing predictions against experimental measurements
4. Using optimisation (Bayesian or gradient-based) to iteratively improve the
   parameter estimates

The included case models a **single laser track on SS316L** at three power
levels (140 W, 200 W, 260 W) and calibrates the simulation to match
the experimentally measured meltpool cross-section geometry.

## Directory Structure

```
AutoCalibrateParameter/
├── main.py                          # Entry point — run this
├── config.yaml                      # Calibration settings (editable)
├── test_environment.py              # Environment check — run this first
├── SingleTrackExperimentalData.csv  # Target experimental data
├── Allrun / Allclean                # Standard OpenFOAM case scripts
├── initial/                         # Initial boundary conditions
├── constant/                        # Material properties & laser config
│   ├── transportProperties          #   SS316L thermo-physical properties
│   ├── LaserProperties              #   Laser model parameters
│   ├── location                     #   Powder bed particle positions
│   ├── timeVsLaserPower             #   Laser power schedule
│   └── timeVsLaserPosition          #   Laser scan path
└── system/                          # Solver settings
    ├── controlDict                  #   Time stepping & output control
    ├── blockMeshDict                #   Mesh definition (300×800×300 µm)
    ├── decomposeParDict             #   Parallel decomposition
    └── ...
```

## Prerequisites

| Requirement | Purpose |
|---|---|
| **OpenFOAM** (v2506 recommended) | CFD solver runtime |
| **laserbeamFoam** solver | Compiled and available in `$PATH` |
| **Python ≥ 3.9** | Calibration framework |
| **ParaView / pvpython** | VTK → CSV meltpool extraction |

### Python Packages

```bash
pip install numpy pandas scipy matplotlib pyyaml joblib scikit-optimize
```

> `scikit-optimize` is only needed for `--method bayes` (the default).
> For `--method gradient`, `scipy` alone is sufficient.

## Quick Start

### 1. Verify your environment

```bash
python test_environment.py
```

This checks that all Python packages are installed, the config file parses
correctly, experimental data loads, and the post-processing dependency chain
(numpy → pandas → joblib → matplotlib) works end-to-end. **Run this first**
— it catches common issues before you spend time on a long calibration run.

### 2. Prepare the OpenFOAM case

```bash
# Source your OpenFOAM environment, e.g.:
of2506            # or: source $WM_PROJECT_DIR/etc/bashrc

# Generate mesh and initial fields
./Allrun
```

### 3. Run the calibration

```bash
# Default: Bayesian optimisation (5 initial points + 5 batches for demo)
python main.py

# Or use gradient-based optimisation
python main.py --method gradient

# See all options
python main.py --help
```

Results will be saved under `./runs/`, with one sub-directory per evaluation
containing the full OpenFOAM case and extracted meltpool geometry.

### 4. Production runs

The default `config.yaml` is tuned for a **quick demonstration** (5 initial
points, 5 batches). For production-quality calibration, edit `config.yaml`:

```yaml
n_initial_points: 20    # more initial exploration
n_batches: 30            # more optimisation iterations
n_proc: 24               # match your CPU count
```

## Configuration Reference

All calibration parameters are set in `config.yaml`:

| Key | Default | Description |
|---|---|---|
| `n_initial_points` | 5 | Initial random samples (Bayesian) |
| `n_batches` | 5 | Optimisation iterations |
| `n_proc` | 4 | MPI processes for OpenFOAM |
| `active_params` | sigma, marangoni, … | Parameters to optimise |
| `fixed_values` | substrate_temp: 300 | Held-constant parameters |
| `output_weights` | [1, 1, 0] | Weights for [width, depth, area] |
| `postproc_script` | (auto-resolved) | Path to meltpool extraction script |

## Experimental Data Format

`SingleTrackExperimentalData.csv` must contain columns:

```
power_W,depth_um,width_um
140,47.53,89.61
200,56.87,111.52
260,63.15,133.1
```

Each row represents one laser power condition with its measured meltpool
depth and width in micrometres.

## Troubleshooting

| Symptom | Check |
|---|---|
| `ImportError: No module named 'skopt'` | `pip install scikit-optimize` |
| `pvpython: command not found` | Install ParaView; add to `$PATH` |
| `meltpool.csv` is empty | Ensure simulation ran to completion |
| Post-processing returns NaN | Check `y_begin_track` / `y_end_track` in config are inside the domain |
| `foam_runner` error | Set `foam_runner: null` if OpenFOAM is already sourced |

## References

- [LaserbeamFoam repository](https://github.com/laserbeamfoam/LaserbeamFoam)
- Bayesian optimisation: [scikit-optimize documentation](https://scikit-optimize.github.io/)
