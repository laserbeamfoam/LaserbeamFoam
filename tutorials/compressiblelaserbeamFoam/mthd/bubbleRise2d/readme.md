# Case Description
Argon bubble rise through liquid metal GaInSn (Galinstan) in the
presence of a horizontally applied magnetic field.

**There are two modes of operation for this case, described in the subsections
below.**

A quasi-incompressible treatment is used presently for GaInSn, while Argon is
modelled as an ideal gas.

This case is based on
_Single bubble rise in GaInSn in a horizontal magnetic field_
(2018) by Richter et al. [1].

## 2.5D Case ($\vec{B} = B \, \vec{e}_z$, page-normal, Richter et al. like)
```
# Enable using:
./Allrun -r -...
```
This case intends to closely match the Richter et al. setup [1], but in a
pseudo-2D sense, as seen in the work of T. F. Flint, M. C. Smith, and
P. Shanthraj [2] where the Richter case was simulated in 2D.

Key points:
- Only one cell in $z$ initially. AMR will refine in the $z$ axis also!
- The $z$-normal walls are all `zeroGradient` for all properties.
This permits the use of a page-normal magnetic field $B_z$ while keeping the
flow predominantly two dimensional.
- Note, in the experiment the channel width in $z$ is only $12$ $\mathrm{mm}$,
  while the bubble equivalent diameter is $5.22$ $\mathrm{mm}$ [1].
  So this case assumes the $z$ channel effects are negligble!

This case may have a slight performance penalty due to the solution of the 3rd
dimension components.

## 2D case ($\vec{B} = B \, \vec{e}_x$, page-horizontal)
```
# Simply avoid the `-r` flag
./Allrun ...
```
This case runs by default. It diverges more signficantly from the Richter et al.
experiment, as it uses a horizontal magnetic field.

This exists for run-time performance by solving for only 2 spatial
dimensions.

# Possible Future Features
**For GaInSn in the Richter et al. case, the temperature variation
is expected to be low, so use of the quasi-incompressible treatment is
reasonable.** That said, for cases when the temperature variation is
significant, please see the comments below.

The GaInSn material should be modelled as a temperature varying,
pressure-quasi-incompressible thermodynamic model. Such as,
$$
    \rho(p,T)
  =
    p/\widetilde{R}_\mathrm{pseudo}
  +
    f(T)
$$
where
- $p/\widetilde{R}$ is a quasi-incompressible term to prevent numerical
  stiffness; and $\widetilde{R}$ is set suitably high, e.g., $10^6$.
- $f(T)$ would be the temperature varying properties.

For example this $f(T)$ could be a lookup table, as seen in Plevachuk et al.
(2014) table 3 [3].

# References
[1] T. Richter et al.,
“Single bubble rise in GaInSn in a horizontal magnetic field,”
International Journal of Multiphase Flow, vol. 104, pp. 32–41, Jul. 2018,
doi: https://doi.org/10.1016/j.ijmultiphaseflow.2018.03.012.

[2] T. F. Flint, M. C. Smith, and P. Shanthraj,
“Magneto-hydrodynamics of multi-phase flows in heterogeneous systems with large
property gradients,” Scientific Reports, vol. 11, no. 1, Sep. 2021,
doi: https://doi.org/10.1038/s41598-021-97177-8.

[3] Y. Plevachuk, V. Sklyarchuk, S. Eckert, Günter Gerbeth, and Rada Novaković,
“Thermophysical Properties of the Liquid Ga–In–Sn Eutectic Alloy,”
Journal of Chemical & Engineering Data, vol. 59, no. 3, pp. 757–763, Feb. 2014,
doi: https://doi.org/10.1021/je400882q.