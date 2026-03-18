# Copilot Instructions

Project context
- Solver split: laserbeamFoam is baseline (Darcy), laserStressFoam includes solid elasticity coupling.
- Linear elasticity is implemented in src/solidElasticity and linked only by laserStressFoam.

Implementation rules
- Keep coupling implicit only; do not add explicit mode.
- Use fSolid = clamp(alpha.metal,0..1) * (1 - clamp(epsilon1,0..1)).
- Use sigma field name sigmaElast (sigma is reserved for surface tension).
- Solve elasticity inside the PIMPLE loop with its own inner iterations.
- [MAYBE?] Add an anchor for displacement to remove rigid-body modes (internal reference or penalty).

Dictionary rules
- Use constant/solidProperties with solidElasticity sub-dictionary:
  - E, nu, alphaT, Tref
  - maxIter, tolerance, relTol
  - solidDamping, fSolidMin
  - coupling = implicit

Tutorials
- laserStressFoam TestStress is the default elasticity test case.
- Add D field in initial/ with zeroGradient on most patches and fixedValue where needed.
- Ensure fvSolution contains D and DFinal solvers.
- Ensure fvSchemes includes div(sigmaElast) and laplacian(DD,D).

Build
- solidElasticity is a lib built via src/Allwmake; laserStressFoam links -lsolidElasticity and includes src/solidElasticity paths.
