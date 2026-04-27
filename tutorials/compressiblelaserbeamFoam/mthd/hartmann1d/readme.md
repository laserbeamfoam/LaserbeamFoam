# Hartmann 1D
**This case is intended as a simple smoke test for the MHD model solvers.**

A planar Hartmann problem solved in steady state. Using a pseudo-incompressible
density model:
$$
\begin{gather*}
    \rho(p,T) = \rho_\mathrm{actual} + \frac{p}{R_\mathrm{pseudo} \, T}
\\
    R_{\mathrm{psuedo}} = 10^{5}
\end{gather*}
$$
(The pseudo-incompressible model is purely to prevent numerical stiffness
when solving for the pressure `p_rgh`).

## Description
This is a simple planar Hartmann flow problem and is similar to a planar
Hagen–Poiseuille flow.

The problem is steady state with 2 spatial dimensions $(x,y)$, with the axial
$x$ direction being fully developed. The $y$ direction is the transverse
direction.

An axial flow, $U_x$, of a MHD fluid is induced by a fixed pressure inlet and
outlet, across a single mesh cell width.
With the inlet and outlet being a `zeroGradient` wall for the velocity
$\vec{U}$.

A magnetic field $H_y$ is applied in the $y$ direction via the top and bottom walls.
Due to the electromagnetic properties of the fluid, namely the magnetic
permittivity $\mu_M$ and the electrical conductivity $\sigma$, an axial magnetic
field $H_x$ can be induced (for non-zero Hartmann numbers).

The Hartmann number may be expressed as
$$
\mathrm{Ha} = B \, L \, \sqrt{ \frac{\sigma}{\mu_f} }
$$
The Hartmann number effectively governs to the resulting solution. In this
case, we set $L, \mu_f = 1$ to allow the user to control the flow by modifying
either $\sigma$ or $B$.
For more details on the Hartmann problem, please see
_Ch. IV: II.2 Hartmann flow_ in [1].

## Usage (Changing Parameters)
- The set magnetic field value $B$ may be changed in `./initial/initialValues`.
- The electrical conductivity may be changed in `./constant/propertiesReference`.
- Due to the usage of a steady state solver, under relaxation factors are
  recommended.

## References
[1] R. Moreau, Magnetohydrodynamics. Springer Nature (Netherlands), 1990.
doi:
[https://doi.org/10.1007/978-94-015-7883-7](https://doi.org/10.1007/978-94-015-7883-7).