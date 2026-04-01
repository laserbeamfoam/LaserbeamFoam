---
name: laserbeamfoam
description: >
  Use this skill for any task involving the LaserbeamFoam repository: editing
  the laserbeamFoam or compressibleLaserbeamFoam solvers, modifying the custom
  ray-tracing or geometric VOF libraries, updating transport or turbulence
  models, adjusting OpenFOAM build files, or changing tutorial cases. Trigger
  whenever the user mentions laserbeamFoam, LaserbeamFoam, laserHeatSource,
  compressibleLaserbeamFoam, laserMeltFoam, geometricVoF, or requests C++ or
  OpenFOAM changes in this codebase.
---

# Codex Guidelines for LaserbeamFoam

This document defines how automated coding changes should be made in this
repository.

## 1) Core Principles

- Make the smallest correct change.
- Preserve existing OpenFOAM-style architecture and naming patterns.
- Prefer consistency with nearby solver or library code over introducing a new
  style.
- Avoid broad refactors unless explicitly requested.
- Keep supported OpenFOAM version compatibility in mind.

## 2) Repository Scope

This repository is a suite of OpenFOAM-based solvers and libraries for
laser-material interaction problems, including:

- `laserbeamFoam`
- `compressibleLaserbeamFoam`
- `laserMeltFoam`
- custom libraries under `src/` for:
  - geometric VOF / interface reconstruction
  - laser ray-tracing heat source modelling
  - transport models
  - turbulence models

The build entry point is the root [`Allwmake`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/Allwmake),
which currently checks for OpenFOAM versions `v2412`, `v2506`, and `v2512`.

## 3) Coding Style Rules

### C++ style

- Follow existing OpenFOAM / repository style in surrounding files.
- Use the same indentation, brace style, comment style, and naming conventions
  as nearby code.
- Keep expressions readable; avoid compact or clever rewrites.
- Prefer explicit local changes over abstraction-heavy redesigns.
- Add comments only when behaviour is non-obvious.

### Headers and includes

- Preserve the existing license/header block format in C++ files.
- Keep include ordering aligned with nearby files.
- Do not change copyright or license text unless explicitly requested.

### Scripts and docs

- Match existing shell style in `Allwmake`, `Allwclean`, `Allrun`, `Allclean`,
  and tutorial helper scripts.
- Keep Markdown concise, practical, and repository-specific.
- Do not rewrite documentation just to restyle it.

## 4) OpenFOAM Conventions to Follow

- Respect runtime type and selection table conventions:
  - `TypeName("...")` in headers
  - `defineTypeNameAndDebug(...)` where appropriate
  - `addToRunTimeSelectionTable(...)` in source files for selectable models
- Preserve dictionary-driven behaviour and runtime configurability.
- Keep field naming and dimension handling consistent with existing solver code.
- Avoid introducing dependencies that break `wmake` / `Allwmake` workflows.
- When adding source files, update the relevant `Make/files` in the same
  module.

Examples of runtime-selectable code in this repository include:

- viscosity models in
  [`src/transportModels/incompressible/viscosityModels`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/src/transportModels/incompressible/viscosityModels)
- transport model factories in
  [`src/transportModels/incompressible`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/src/transportModels/incompressible)
- turbulence model factories in
  [`src/turbulenceModel`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/src/turbulenceModel)

Rules when adding a new runtime-selectable class:

- Add `TypeName("...")` in the header.
- Register the type in the matching source file when the module uses a runtime
  selection table.
- Add the source file to the relevant `Make/files`.
- Ensure the dictionary `type` string matches the registered class name.

## 5) Solver and Library Patterns

When editing the main solvers, preserve the existing split-header structure:

- solver entry point in `*.C`
- field creation in `createFields.H`
- momentum equation in `UEqn.H`
- pressure equation in `pEqn.H`
- energy equation in `TEqn.H`
- property updates in files such as `updateProps.H`, `updateKappaCp.H`, or
  `update.H`
- VOF controls under solver-local `MULES/`, `isoAdvector/`, or `VoF/`
  directories

When editing laser deposition logic:

- Prefer keeping ray-tracing behaviour in the `laserHeatSource` library rather
  than pushing detail into solver loops.
- Preserve dictionary-based laser configuration from `constant/LaserProperties`
  and time-series inputs such as `timeVsLaserPosition` and
  `timeVsLaserPower`.
- Be careful with changes that affect ray-path VTK output or deposition field
  coupling into `TEqn.H`.

When editing geometric VOF code:

- Follow existing separation between reconstruction schemes, sampled
  interfaces, surface iterators, and advection code.
- Avoid mixing solver-specific logic into the generic `src/geometricVoF`
  library.

## 6) Build and Test Workflow

Preferred workflow:

1. Read nearby code and match local patterns.
2. Make the minimum patch necessary.
3. Update `Make/files` or `Make/options` if required.
4. Build the narrowest relevant target first when possible.
5. Run the most relevant tutorial, utility, or compile check available.

Typical build entry points:

- root build: [`Allwmake`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/Allwmake)
- library build: [`src/Allwmake`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/src/Allwmake)
- application build: [`applications/Allwmake`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/applications/Allwmake)
- compressible solver-specific build:
  [`applications/solvers/compressibleLaserbeamFoam/Allwmake`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/applications/solvers/compressibleLaserbeamFoam/Allwmake)

Relevant regression surfaces usually include tutorial cases under
[`tutorials`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/tutorials),
especially:

- `tutorials/laserbeamFoam/...`
- `tutorials/compressiblelaserbeamFoam/...`
- `tutorials/multiComponentlaserbeamFoam/...`
- `tutorials/laserMeltFoam/...`

## 7) Minimise Changes

- Only modify files necessary for the requested task.
- Do not reformat unrelated code.
- Do not rename files or symbols unless required.
- Do not change solver behaviour outside requested scope.
- Prefer targeted fixes over opportunistic cleanup.

Before finalizing, verify:

- build impact is localized
- only relevant files changed
- no solver branch or library registration was accidentally broken
- dictionary interfaces remain compatible unless the user asked to change them

## 8) What to Avoid

- Large-scale refactors without explicit request.
- Moving code between solver and library layers without clear need.
- Introducing new coding conventions inconsistent with the repo.
- Silent changes to case dictionary interfaces.
- Updating tutorial inputs as a side effect of solver edits unless necessary.

## 9) Reference Files

- Main incompressible solver:
  [`applications/solvers/laserbeamFoam/laserbeamFoam.C`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/applications/solvers/laserbeamFoam/laserbeamFoam.C)
- Main compressible solver:
  [`applications/solvers/compressibleLaserbeamFoam/compressibleLaserbeamFoam.C`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/applications/solvers/compressibleLaserbeamFoam/compressibleLaserbeamFoam.C)
- Laser heat source library:
  [`src/laserHeatSource/laserHeatSource.H`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/src/laserHeatSource/laserHeatSource.H)
  [`src/laserHeatSource/laserHeatSource.C`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/src/laserHeatSource/laserHeatSource.C)
- Ray particle implementation:
  [`src/laserHeatSource/laserRayParticle.H`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/src/laserHeatSource/laserRayParticle.H)
  [`src/laserHeatSource/laserRayParticle.C`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/src/laserHeatSource/laserRayParticle.C)
- Geometric VOF core:
  [`src/geometricVoF/Make/files`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/src/geometricVoF/Make/files)
  [`src/geometricVoF/reconstructionSchemes/reconstructionSchemes.H`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/src/geometricVoF/reconstructionSchemes/reconstructionSchemes.H)
- Incompressible transport models:
  [`src/transportModels/incompressible/Make/files`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/src/transportModels/incompressible/Make/files)
- Turbulence model entry points:
  [`src/turbulenceModel/turbulenceModelNew.H`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/src/turbulenceModel/turbulenceModelNew.H)
  [`src/turbulenceModel/turbulenceModelNew.C`](/Users/philipc/code/laserbeamfoamProjects/LaserbeamFoam/src/turbulenceModel/turbulenceModelNew.C)
