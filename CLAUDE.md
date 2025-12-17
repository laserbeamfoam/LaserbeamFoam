# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

LaserbeamFoam is an OpenFOAM-based suite of solvers for simulating high-energy-density laser-substrate interactions in advanced manufacturing processes (laser welding, drilling, powder bed fusion, selective laser melting). The repository contains three main solvers:

- **laserbeamFoam**: VOF solver for two-phase incompressible flow with laser heat source and ray-tracing
- **laserMeltFoam**: Extension with additional liquid fraction handling
- **compressibleLaserbeamFoam**: Multi-component metallic substrates with vapor/condensed phase transitions

## Build System

This is an OpenFOAM extension that uses the OpenFOAM/FOAM wmake build system.

### Requirements
- OpenFOAM v2506 (for OpenFoam_com_main branch) or OpenFOAM10 (for Openfoam_org_main branch)
- Optional: LIGGGHTS® for DEM powder bed simulations

### Build Commands

```bash
# Build everything (libraries + applications)
./Allwmake -j

# Clean everything
./Allwclean

# Build logs are written to log.Allwmake in each component directory
# Check for errors after build
find . -name log.Allwmake | xargs grep " Error "
find . -name log.Allwmake | xargs grep " Stop."
```

The build process:
1. Compiles libraries in `src/` (geometricVoF, laserHeatSource, transportModels, turbulenceModel)
2. Compiles applications in `applications/` (solvers and utilities)

## Code Architecture

### Source Libraries (`src/`)

**geometricVoF**: Geometric Volume-of-Fluid interface tracking
- `isoAdvection`: Iso-surface advection scheme for sharp interface tracking
- `reconstructionSchemes`: Interface reconstruction methods (PLIC, RDF)
- `cellCuts`: Geometric cutting algorithms for interface cells
- `surfaceIterators`: Interface surface iteration utilities

**laserHeatSource**: Ray-tracing laser beam implementation
- `laserHeatSource.H/C`: Main heat source class with Fresnel equations and multi-reflection ray-tracing
- `compactRay.H/C`: Individual ray representation and tracking
- Reads laser properties from `constant/LaserProperties`
- Writes ray paths to `VTK/rays_<LASER_NAME>_<TIME>.vtk` for ParaView visualization

**transportModels**: Two-phase mixture properties with phase-change
- Handles fusion/melting, vaporization, and associated latent heat
- Temperature-dependent properties (surface tension, Marangoni effects)
- Momentum damping due to solidification

**turbulenceModel**: Turbulence modeling extensions for two-phase laser-material interaction

### Solver Structure

All solvers follow a similar PIMPLE algorithm structure:

1. Dynamic time-step adjustment for stability
2. Phase equation sub-cycle (MULES or isoAdvector)
3. Update interface location for heat source
4. Update fluid properties
5. Ray-tracing heat source application at surface
6. PIMPLE loop:
   - Form U equation
   - Energy transport loop (solve T, update phase fractions, latent heat)
   - PISO pressure-velocity correction

Key solver files use `#include` directives for equation components (e.g., `UEqn.H`, `TEqn.H`, `pEqn.H`, `alphaEqnSubCycle.H`) located in solver-specific subdirectories.

### Applications

**Solvers** (`applications/solvers/`):
- `laserbeamFoam/`: Main incompressible two-phase solver
- `laserMeltFoam/`: With liquid fraction tracking
- `compressibleLaserbeamFoam/`: Multi-component with compressibility

**Utilities** (`applications/utilities/`):
- `setSolidFraction`: Initializes metal volume fraction from DEM particle locations

## Running Tutorials

Tutorial cases are in `tutorials/<solver_name>/`. Each has:
- `Allrun`: Setup and run script
- `Allclean`: Clean case
- `initial/`: Initial field templates
- `constant/`: Case properties including `LaserProperties`, material properties
- `system/`: discretization and solver settings

### Standard Tutorial Workflow

```bash
cd tutorials/laserbeamFoam/<case_name>
./Allrun       # Copies initial/ to 0/, runs blockMesh, setFields, solver
./Allclean     # Removes generated files
```

Typical Allrun steps:
```bash
cp -r initial 0
blockMesh                    # Create mesh
setFields                    # Initialize fields
laserbeamFoam               # Run solver (serial)

# For parallel:
# decomposePar
# mpirun -np 6 laserbeamFoam -parallel &> log.laserbeamFoam
```

### LPBF Cases with DEM

LPBF (Laser Powder Bed Fusion) tutorials use DEM-generated powder beds:
1. `DEM_small/Allrun` runs LIGGGHTS to generate `post/location` (particle positions)
2. `setSolidFraction` reads `constant/location` to initialize metal volume fraction field
3. Pre-generated powder beds are included, so LIGGGHTS is optional

## Testing

```bash
# Run all tutorial tests
cd tutorials
./Alltest
```

## Visualization

Ray paths: ParaView can load `VTK/rays_<LASER_NAME>.vtk.series` with correct time synchronization. Use "Tube" filter or increase "Line Width" for better visibility.

## Key Configuration Files

- `constant/LaserProperties`: Laser parameters (power, beam radius, scanning path)
- `constant/transportProperties`: Phase properties, surface tension, Marangoni coefficients
- `constant/thermophysicalProperties`: Thermal properties, melting/boiling points, latent heats
- `system/controlDict`: Time stepping, output settings
- `system/fvSchemes`: Numerical schemes
- `system/fvSolution`: Solver settings and PIMPLE parameters

## Branch Structure

- `OpenFoam_com_main`: Compatible with OpenFOAM v2506 (openfoam.com)
- `Openfoam_org_main`: Compatible with OpenFOAM10 (openfoam.org)
