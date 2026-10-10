# microstructureFoam

## Overview

microstructureFoam is an OpenFOAM-based solver for simulating microstructure evolution during laser-based additive manufacturing processes such as powder bed fusion (PBF) and directed energy deposition (DED). It couples a two-phase incompressible Volume of Fluid (VoF) description of the metallic substrate and shielding gas with a multi-order-parameter phase-field model for grain nucleation, growth, competitive selection, and coarsening.

The solver includes a ray-tracing laser heat source model that computes Fresnel absorption at the metal-gas interface, accounting for multiple reflections within surface depressions such as keyholes. The phase-field model captures the liquid-solid phase transformation driven by the local temperature field, and reproduces key grain evolution phenomena including epitaxial growth from the substrate and powder particles, heterogeneous nucleation in undercooled liquid, grain coarsening in heat-affected zones, re-melting and re-growth in overlapping zones between adjacent tracks and layers, and competitive grain growth with crystallographic anisotropy.

The underlying phase-field formulation follows the model of Yang, Wang and Yan (npj Computational Materials, 2021), which uses non-conserved order parameters to describe both the liquid-solid transformation and grain orientations. Two classes of phase-field variable are introduced: a solid fraction parameter ξ (`xi`) that distinguishes liquid (ξ=0) from solid (ξ=1), and a set of orientation order parameters η_i (`n.0`, `n.1`, ..., `n.N`) that each identify a grain of a specific crystallographic orientation, where η_i=1 inside the i-th grain and 0 elsewhere. The liquid fraction ε₁ (`epsilon1`) is derived as ε₁ = 1 − ξ after each phase-field solve.

The governing equations are time-dependent Ginzburg-Landau equations that minimise the system's total free energy, consisting of a phase free energy (coupling ξ to the local temperature through the liquidus), a grain interaction free energy (penalising grain overlap and coupling grain identity to the solid/liquid state), and gradient energy terms that regularise the liquid-solid and grain boundary interfaces.

The grain boundary mobility is temperature-dependent via an Arrhenius relation, and anisotropic grain boundary energy can be included to promote preferential growth along crystallographic directions (e.g. ⟨001⟩ for cubic metals). Nucleation is handled through a site-saturation model: nucleation sites are distributed randomly throughout the domain at initialisation with activation undercoolings drawn from a Gaussian distribution. Sites activate when the local undercooling exceeds their assigned threshold and the local region is free of existing grains, introducing a new order parameter at that location.


## Installation

The current version of the code is developed for [OpenFOAM v2506](https://www.openfoam.com/news/main-post/openfoam-v2506) (ESI/OpenCFD distribution). The code has been developed and tested using an Ubuntu installation, but should work on any operating system capable of installing OpenFOAM.

To install, first follow the instructions at the [OpenFOAM v2506 installation page](https://develop.openfoam.com/Development/openfoam/-/wikis/precompiled/debian) to install the OpenFOAM v2506 libraries.

Then navigate to a working folder in a shell terminal, clone the git code repository, and build:

```bash
git clone https://github.com/micmog/laserbeamFoam.git laserbeamFoam
cd laserbeamFoam/solver
wclean
wmake
```

The installation can be tested using the tutorial cases described below.

### Building the TesselateFoam utility

An optional pre-processing utility, TesselateFoam, populates the initial grain structure using Voronoi tessellation. It assigns each mesh cell to the nearest randomly generated seed point, setting the corresponding order parameter to 1.0 in that cell. To build it:

```bash
cd laserbeamFoam/utilities/TesselateFoam
wclean
wmake
```


## Tutorial Cases

To run any of the tutorials in serial mode:

```bash
# Clean any old simulation files
rm -r 0* 1* 2* 3* 4* 5* 6* 7* 8* 9*

# Copy initial conditions
cp -r initial 0

# Generate the mesh
blockMesh

# Set initial field values
setFields

# (Optional) Initialise grain structure via Voronoi tessellation
TesselateFoam

# Run the solver
microstructureFoam
```

For parallel deployment using MPI, after `setFields`:

```bash
decomposePar
mpirun -np 6 microstructureFoam -parallel > log &
```

### Powder-Bed Fusion Example

This tutorial demonstrates a laser powder bed fusion simulation with a moving Gaussian heat source scanning across a metallic substrate with a shielding gas atmosphere. The ray-tracing model deposits energy at the metal-gas interface, driving melting, fluid flow (including Marangoni convection and recoil pressure), solidification, and grain evolution. Nucleation sites activate as the melt pool solidifies and undercooling develops. The phase-field equations then govern competitive grain growth behind the solidification front.

### Directional Solidification Example

This tutorial demonstrates solidification under a prescribed temperature gradient without a laser source, allowing the study of grain selection and competitive growth in isolation. It is useful for validating the phase-field model parameters and observing the columnar-to-equiaxed transition (CET) as a function of temperature gradient and cooling rate.


## Description of Case Files

### transportProperties

General transport and phase properties for the two-phase system.

| Parameter Name | Value Type | Description | Unit |
| --- | --- | --- | --- |
| `phases` | List\[Word\] | Names of the two phases (e.g. metal and gas). A separate `physicalProperties.[PHASE]` dictionary is read for each phase | - |
| `sigma` | scalar | Surface tension coefficient between the two phases | N/m |
| `dsigmadT` | scalar | Temperature coefficient of surface tension (Marangoni coefficient, dσ/dT) | N/(m·K) |
| `p0` | scalar | Atmospheric/reference pressure, used in the recoil pressure and evaporation models | Pa |
| `Tvap` | scalar | Vaporisation temperature of the metallic phase | K |
| `Mm` | scalar | Molar mass of the metallic phase | kg/mol |
| `LatentHeatVap` | scalar | Latent heat of vaporisation | J/kg |
| `emL` | scalar | Emissivity of the liquid phase (for radiative heat loss) | - |
| `emS` | scalar | Emissivity of the solid phase (for radiative heat loss) | - |
| `TRef` | dimensionedScalar | Reference temperature for radiative heat transfer | K |
| `interfaceTrackingScheme` | word | Interface capturing method: `MULES` or `isoAdvector` | - |

### physicalProperties.\[PHASE\]

Per-phase material properties. One dictionary is required for each phase listed in `transportProperties`.

| Parameter Name | Value Type | Description | Unit |
| --- | --- | --- | --- |
| `transportModel` | word | Viscosity model (typically `Newtonian`) | - |
| `nu` | scalar | Kinematic viscosity | m²/s |
| `rho` | scalar | Density | kg/m³ |
| `Tsolidus` | scalar | Solidus temperature | K |
| `Tliquidus` | scalar | Liquidus temperature | K |
| `LatentHeat` | scalar | Latent heat of solidification/melting | J/kg |
| `beta` | scalar | Thermal expansion coefficient for buoyancy (Boussinesq approximation) | 1/K |
| `poly_kappa` | Polynomial\<8\> | Polynomial coefficients for temperature-dependent thermal conductivity κ(T) | W/(m·K) |
| `poly_cp` | Polynomial\<8\> | Polynomial coefficients for temperature-dependent specific heat capacity cp(T) | J/(kg·K) |
| `elec_resistivity` | scalar | Electrical resistivity (used in Fresnel absorption calculation via the Drude model) | Ω·m |

### nucleationDict

Controls the heterogeneous nucleation model. Nucleation sites are randomly distributed throughout the mesh at startup, with activation undercoolings sampled from a Gaussian distribution.

| Parameter Name | Value Type | Description | Unit |
| --- | --- | --- | --- |
| `Tu_mean` | dimensionedScalar | Mean of the Gaussian activation undercooling distribution | K |
| `Tu_stdev` | dimensionedScalar | Standard deviation of the Gaussian activation undercooling distribution | K |
| `n_max` | dimensionedScalar | Maximum volumetric density of nucleation sites. For 2D meshes this is automatically converted to an areal density | sites/m³ |
| `nucDistFactor` | scalar | Minimum spacing between nucleation sites, expressed as a multiple of the phase-field interface width max(lp, lg) | - |
| `maxItersNucSet` | scalar | Maximum iterations allowed when placing each nucleation site (to prevent infinite loops if the domain is saturated) | - |
| `stopToCheckNucSites` | Switch | If true, the simulation stops after setting nucleation sites so that the `nucTu` and `nucSite` fields can be inspected in ParaView before running | true/false |

### LaserProperties

Controls the ray-tracing laser heat source. The laser is modelled as a Gaussian beam with rays traced through the domain, depositing energy at the metal-gas interface via Fresnel absorption computed from the Drude model.

| Parameter Name | Value Type | Description | Unit |
| --- | --- | --- | --- |
| `V_incident` | vector | Initial incidence direction vector for the laser beam | - |
| `HS_a` | scalar | Laser beam radius (1/e² definition) | m |
| `HS_bg` | scalar | Laser beam centre offset in the x-direction | m |
| `HS_lg` | scalar | Laser beam centre offset in the z-direction | m |
| `HS_velocity` | scalar | Travel velocity of the laser heat source in the z-direction | m/s |
| `HS_Q` | scalar | Laser power | W |
| `wavelength` | scalar | Wavelength of the incident laser radiation | m |
| `e_num_density` | scalar | Free electron number density of the substrate material (Drude model) | m⁻³ |
| `elec_resistivity` | scalar | Electrical resistivity of the substrate material (Drude model) | Ω·m |
| `PowderSim` | Switch | Enable powder bed simulation mode (uses `epsilon1mask` for Darcy damping to prevent spurious currents smearing out partially melted particles) | true/false |
| `useLocalSearch` | Switch | Use local cell-to-cell search for ray tracing (faster than global N² search) | true/false |
| `maxLocalSearch` | label | Maximum number of cells to check during a local search before falling back to global search | - |
| `debug` | Switch | Enable verbose ray-tracing debug output | true/false |

### PhaseFieldProperties

Controls the multi-order-parameter phase-field model for grain evolution.

| Parameter Name | Value Type | Description | Unit |
| --- | --- | --- | --- |
| `N_Seeds` | label | Number of initial grain seeds (set by TesselateFoam or read from existing n.i fields). The total number of order parameters is N_Seeds + number of nucleation sites | - |
| `gamma_pf` | scalar | Phase-field model parameter γ controlling the grain boundary energy relative to bulk grain interaction. Determines the energy penalty for grain overlap | - |
| `sigma_P` | scalar | Liquid-solid interfacial energy σ_p | J/m² |
| `sigma_G` | scalar | Grain boundary energy σ_g | J/m² |
| `lp` | scalar | Liquid-solid interface diffuse width l_p | m |
| `lg` | scalar | Grain boundary diffuse width l_g | m |
| `deltaf_p` | scalar | Maximum height of the free energy barrier at the liquid-solid interface Δf_p | - |
| `deltaf_g` | scalar | Maximum height of the free energy barrier at grain boundaries Δf_g | - |
| `D0` | scalar | Pre-exponential coefficient for grain boundary mobility D_g = D_0 exp(-Q_g / RT) | m³/(J·s) |
| `Qg` | scalar | Activation energy for grain boundary mobility | J/mol |
| `Aniso_strength` | scalar | Strength of crystallographic anisotropy in grain boundary energy (ε₀ in the anisotropy function). Set to 0 to disable anisotropy | - |
| `grainNumThreshold` | scalar | Minimum value of an order parameter η_i for a cell to be considered as belonging to grain i (used for grain number visualisation) | - |
| `Xmin`, `Xmax` | scalar | Bounding box limits in x for Voronoi tessellation (TesselateFoam) | m |
| `Ymin`, `Ymax` | scalar | Bounding box limits in y for Voronoi tessellation (TesselateFoam) | m |
| `Zmin`, `Zmax` | scalar | Bounding box limits in z for Voronoi tessellation (TesselateFoam) | m |
| `write_ni_all` | Switch | Write all order parameter fields n.i at write times (including inactive/empty fields). Default: true | true/false |
| `write_ni_active` | Switch | Write only active (non-zero) order parameter fields. Only used if `write_ni_all` is false. Default: true | true/false |

### g

| Parameter Name | Value Type | Description | Unit |
| --- | --- | --- | --- |
| `dimensions` | dimensionSet | Dimension set for gravitational acceleration | m/s² |
| `value` | vector | Direction and magnitude of gravitational acceleration | m/s² |

### fvSolution / MELTING

| Parameter Name | Value Type | Description | Unit |
| --- | --- | --- | --- |
| `damperSwitch` | bool | Enable density-ratio damping for surface force and thermal source terms. Helps stability with large density ratios between metal and gas | true/false |


## Algorithm

The solver initialises the mesh, reads fields and boundary conditions, distributes nucleation sites, initialises the grain structure (optionally via Voronoi tessellation), and selects the turbulence model. The main time loop then proceeds as follows:

1. **Update Δt** — The time step is dynamically adjusted to satisfy Courant number constraints for both the bulk flow and the interface advection.

2. **PIMPLE pressure-velocity loop** — The following steps are repeated for a user-specified number of outer correctors:

    a. **VOF equation sub-cycle** — The phase-fraction equation for α₁ (metal volume fraction) is solved for a user-defined number of sub-time-steps using either MULES or isoAdvector (selected via `interfaceTrackingScheme`), maintaining α₁ ∈ [0, 1].

    b. **Update material properties** (`updateProps.H`) — Mixture properties (cp, κ, β, latent heat, solidus/liquidus temperatures) are updated based on the current α₁ and temperature fields via temperature-dependent polynomials. The filtered interface normal n̂, liquid fraction mask ε₁, and Boussinesq buoyancy factor are recomputed. If anisotropy is enabled (`Aniso_strength > 0`), the orientation-dependent grain boundary mobility is calculated from the order parameter gradients rotated into each grain's crystallographic frame using stored quaternions.

    c. **Laser heat source** — Ray tracing computes the volumetric energy deposition at the metal-gas interface. Rays originate from the beam footprint on the top boundary, propagate through the domain along the incidence direction, and deposit energy via Fresnel absorption (computed from the Drude model using the local surface normal and material optical properties). Reflected rays continue propagating until their energy is negligible, capturing multiple reflections within surface depressions.

    d. **Momentum equation** (`UEqn.H`) — Solves for velocity including Darcy damping in the mushy zone (penalising flow in solid regions via the Kozeny-Carman relation), Marangoni surface tension gradient forces (zeroed below T_liquidus), recoil pressure from vaporisation, and buoyancy.

    e. **Energy equation** (`TEqn.H`) — Solves the enthalpy-based energy equation in a single pass. Source terms include latent heat release/absorption from solidification/melting (via ∂ε₁/∂t), evaporative cooling, linearised radiative heat loss (Stefan-Boltzmann), and the laser deposition field.

    f. **Pressure correction** (`pEqn.H`) — PISO-style pressure-velocity coupling corrects the pressure and velocity fields to satisfy continuity, with non-orthogonal corrections.

    g. **Nucleation** (`nucleate.H`) — Active nucleation sites in liquid metal regions are checked against the local undercooling (T_liquidus − T). If the undercooling exceeds the site's activation threshold and no existing grain occupies the cell (Ση_i ≈ 0), a new order parameter is set to 1.0 at that location and ξ is set to 1.0, representing a newly nucleated grain. Nucleation is synchronised across MPI processors by cycling through ranks sequentially and broadcasting the latest field index after each processor nucleates.

    h. **Phase-field equations** (`PFEqns.H`) — The Ginzburg-Landau equations are solved for all active order parameters η_i. The variational derivative S_i of the grain interaction free energy is computed using an O(N) algorithm that precomputes Σ η_j² to avoid the naive O(N²) pairwise summation. The solid fraction ξ is then solved, coupling to the temperature field through a hyperbolic tangent function of T/T_liquidus. After solving, ξ is clamped to [0, 1].

    i. **Update liquid fraction** — After the phase-field solve, the liquid fraction is updated as ε₁ = max(min(1 − ξ, 1), 0).

3. **Post-PIMPLE updates** — The melt history field is incremented for cells that are simultaneously metallic (α₁ > 0.5) and liquid (ε₁ > 0.5).

4. **Write fields** — At write times, the grain number field (`grainNum`) is recomputed by finding the dominant order parameter in each cell, the crystallographic orientation quaternions (`qw`, `qv`) are populated from the stored rotation quaternions for each grain, and the order parameter fields are written according to the `write_ni_all`/`write_ni_active` settings. Ray path VTK files are also written for visualisation.


## Key Fields

| Field | Description |
| --- | --- |
| `alpha.metal` (α₁) | Metal volume fraction (1 = metal, 0 = gas) |
| `T` | Temperature |
| `U` | Velocity |
| `p_rgh` | Pressure minus hydrostatic component |
| `xi` (ξ) | Solid fraction from the phase-field model (1 = solid, 0 = liquid). Primary variable solved by PFEqns |
| `epsilon1` (ε₁) | Liquid fraction, derived as 1 − ξ after each phase-field solve |
| `epsilon1mask` | Filtered liquid fraction used for Darcy damping (zeroed where the cell-averaged ε₁ ≤ 0.95) |
| `n.0, n.1, ..., n.N` (η_i) | Phase-field order parameters identifying individual grains |
| `grainNum` | Integer field indicating which order parameter dominates each cell (-1 = no grain) |
| `qw`, `qv` | Quaternion components representing the crystallographic orientation of each grain |
| `Deposition` | Laser energy deposited per unit volume at the metal-gas interface |
| `Ray_Q` | Volumetric laser energy from ray tracing |
| `Ray_Number` | Ray identifier field (for debugging ray paths) |
| `DC` | Darcy coefficient penalising flow in the mushy/solid zone |
| `Marangoni` | Marangoni (thermocapillary) force field |
| `pVap` | Recoil pressure from vaporisation |
| `Qv` | Evaporative cooling volumetric source |
| `nucSite` | Nucleation site status (-1 = not a site, 0 = inactive site, 1 = active site) |
| `nucTu` | Activation undercooling assigned to each nucleation site |
| `GBs` | Grain boundary indicator field (Σ η_i²) |
| `Aniso_Mob` | Anisotropic mobility field (weighted average of orientation-dependent mobility) |
| `meltHistory` | Cumulative count of time steps each cell has been simultaneously metallic and liquid |
| `gradT` | Temperature gradient |
| `rhok` | Boussinesq buoyancy factor (1 − ε₁·β·(T − T_solidus)) |


## License

OpenFOAM, and by extension the microstructureFoam solver, is licensed free and open source only under the [GNU General Public Licence version 3](https://www.gnu.org/licenses/gpl-3.0.en.html). One reason for OpenFOAM's popularity is that its users are granted the freedom to modify and redistribute the software and have a right of continued free use, within the terms of the GPL.


## Acknowledgements


## Citing This Work

If you use microstructureFoam in your research, please consider citing the following:

ADD microstructureFoam SoftwareX paper here!

The phase-field model implementation is based on:

> M. Yang, L. Wang, W. Yan, "Phase-field modeling of grain evolutions in additive manufacturing from nucleation, growth, to coarsening", *npj Computational Materials*, **7**, 56 (2021).

Coupled with the thermal-fluid-dynamics solver in laserbeamFoam!


## References

1. M. Yang, L. Wang, W. Yan, "Phase-field modeling of grain evolutions in additive manufacturing from nucleation, growth, to coarsening", *npj Computational Materials*, **7**, 56 (2021).
2. N. Moelans, B. Blanpain, P. Wollants, "Quantitative analysis of grain boundary properties in a generalized phase field model for grain growth in anisotropic systems", *Physical Review B*, **78**, 024113 (2008).
3. G. P. Krielaart, S. V. D. Zwaag, "Simulations of pro-eutectoid ferrite formation using a mixed control growth model", *Materials Science and Engineering A*, **246**, 104–116 (1998).
4. T. DebRoy et al., "Additive manufacturing of metallic components – process, structure and properties", *Progress in Materials Science*, **92**, 112–224 (2018).