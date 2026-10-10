CorrectPhi
(
    U,
    phi,
    p_rgh,
    surfaceScalarField("rAUf", fvc::interpolate(rAU())),
    vDot, // geometricZeroField(),
    pimple
);

Info<<"Corrected phi for Dilation"<<endl;


#include "continuityErrs.H"
