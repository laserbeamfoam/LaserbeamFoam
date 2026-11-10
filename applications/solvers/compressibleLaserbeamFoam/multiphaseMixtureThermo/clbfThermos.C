/*---------------------------------------------------------------------------*\
    clbfThermo: temperature-dependent viscosity and conductivity for the
    condensed phases of compressibleLaserbeamFoam
\*---------------------------------------------------------------------------*/

#include "rhoThermo.H"
#include "makeThermo.H"

#include "specie.H"
#include "adiabaticPerfectFluid.H"
#include "hConstThermo.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"

#include "tabulatedTransport.H"
#include "polynomialTransport.H"

#include "heRhoThermo.H"
#include "pureMixture.H"

namespace Foam
{

// ESI pairs adiabaticPerfectFluid only with constTransport (rhoThermos.C), so mu(T) and kappa(T) need these.
// The energy (hConst) and the equation of state stay as in the cases; only the transport becomes a function of T
makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    tabulatedTransport,
    sensibleInternalEnergy,
    hConstThermo,
    adiabaticPerfectFluid,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    polynomialTransport,
    sensibleInternalEnergy,
    hConstThermo,
    adiabaticPerfectFluid,
    specie
);

} // End namespace Foam

// ************************************************************************* //
