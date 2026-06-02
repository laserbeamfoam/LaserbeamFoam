## Comments to LaserbeamFoam devs:
Todo:
- Implement a temperature varying, pressure-quasi-incompressible
  thermodynamic model? Richter et al. and their citiations report
  GaInSn to have temperature dependent properties, e.g., $\phi(T)$.

Actual read me:
# Case Description
Argon bubble rise through liquid metal GaInSn (Galinstan) in the
presence of a horizontally applied magnetic field. The case is
setup in 2-spatial dimensions $(x,y)$.

A quasi-incompressible treatment is used presently.

This case is based on
_Single bubble rise in GaInSn in a horizontal magnetic field_
(2018) by Richter et al.
(https://doi.org/10.1016/j.ijmultiphaseflow.2018.03.012). However,
for this 2D case, we deviated by apply the magnetic field in the $x$
direction, and not the $z$ direction.
