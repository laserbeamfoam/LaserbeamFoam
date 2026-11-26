// Flux correction after a mesh change (and at a fresh start). divUCorrectPhi, set by the including
// file, is the divergence the corrected flux must have
CorrectPhi
(
    U,
    phi,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU())),
    divUCorrectPhi,
    pimple
);

#include "continuityErrs.H"
