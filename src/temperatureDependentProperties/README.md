# temperatureDependentProperties

Temperature-dependent cp, k, viscosity, electrical conductivity and magnetic
permeability for laserbeamFoam. Each property is any `Function1` of T in the
phase sub-dictionary of `transportProperties`, sampled once at start-up onto a
uniform T grid (`temperatureTable`, default 0–10000 K every 1 K, 80 kB) and
evaluated per cell with one interpolation (~2 ns, the same as the old
`Polynomial<8>`; exact for linear polynomials and tables with nodes on the
grid, otherwise `deltaT²|f''|/8`).

```
propertyTables { Tmin 0; Tmax 10000; deltaT 1; }    // optional

metal
{
    transportModel  temperatureTable;                      // nu(T)
    temperatureTableCoeffs { nu table ((1928 9.9e-7) (3500 3.4e-7)); }   // m2/s
    rho               4430;
    cp                table ((300 546) (1928 831) (5000 831));            // J/kg/K
    kappa             polynomial ((3.4 0) (0.0125 1));                    // W/m/K
    Elec_conductivity table ((300 5.9e5) (1928 5.0e5) (5000 4.4e5));      // S/m
    Mag_permitivity   1.256e-6;                                           // only with electromagnetics on
    ...
}
gas
{
    transportModel  Newtonian;  nu 1.4e-5;
    cp      520;
    kappa   tableFile; kappaCoeffs { file "constant/kArgon.dat"; }
    ...
}
```

Fallbacks so existing cases run unchanged: `poly_cp` / `poly_kappa`
(`Polynomial<8>`, a0 + a1 T + …) when `cp` / `kappa` are absent; the
deprecated constant `elec_resistivity` when `Elec_conductivity` is absent.
A bare number is a constant Function1. Tables narrower than the cache are
extended by the Function1 `outOfBounds` (default clamp).

## Where each property goes

| property | used by | evaluated |
|---|---|---|
| cp, kappa | energy equation (`updateKappaCp.H`), thermalDamper | once per outer corrector at the previous T, as before; cp mass-weighted, kappa volume-weighted |
| nu | `mixture.correct()` → momentum equation | once per outer corrector; Darcy damping still owns the solid |
| Elec_conductivity (metal) | laser absorptivity: `electrical_resistivity = 1/sigma_metal(T)` everywhere (interface cells must look like metal, no gas blending) | same loop |
| Elec_conductivity, Mag_permitivity (both phases) | MHD, when `electromagnetics` is on: alpha-blended fields handed to `mthdModel::setTransportProperties` | same loop |

The H equation needs no change for σ(T): −∇·(η∇H) + ∇·(η(∇H)ᵀ) is
identically ∇×(η∇×H) for a spatially varying η = 1/σ, and the Joule source
is already (1/σ)|J|². `mthdModel` only gains the setter; the old
dictionary/accumulate interface used by compressibleLaserbeamFoam is intact.

The energy equation is already in the advective form, so a variable cp
enters exactly as ρ cp(T) DT/Dt. One deliberate fix rides along: `cp` and
`kappa` were `calculated` fields whose patch values stayed at 0, which made
every boundary adiabatic and enthalpy-flux-free whatever the T condition
said; they are now `zeroGradient`.

The solver prints every table at start-up ("Temperature dependent
properties of metal: cp: Function1 table, cached … values: 300K:546 …").

## Files

```
src/temperatureDependentProperties/   new library (sibling of src/MTHD, built by src/Allwmake)
src/MTHD/mthdModel.{H,C}              + setTransportProperties(), electricalConductivity()
src/Allwmake, src/Allwclean           + one line each
applications/solvers/laserbeamFoam/
    createFields.H   tables, electrical fields, mthd created earlier, zeroGradient cp/kappa
    updateKappaCp.H  table evaluation, resistivity and MHD fields
    updateProps.H    mthd->setTransportProperties after updateKappaCp.H
    TEqn.H           thermalDamper uses the tables
    laserbeamFoam.C  #include "temperatureTable.H"
    Make/options     include path + -ltemperatureDependentProperties
```
compressibleLaserbeamFoam is untouched; when its turn comes, the thermo
side goes through its rhoThermo (`hPolynomial`/`hTabulated`,
`polynomial`/`tabulated` transport) and the MHD side through
`mthdModel::accumulateTransportProperties` with per-phase fields.
