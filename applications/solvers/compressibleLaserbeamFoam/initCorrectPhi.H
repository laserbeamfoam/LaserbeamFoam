tmp<volScalarField> rAU;

if (correctPhi)
{
    rAU = new volScalarField
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rAU", dimTime/dimDensity, 1)
    );

    if (!bRestartFirstLoop)
    {
        #include "correctPhi.esi.H"
    }
}
else if (!bRestartFirstLoop)
{
    CorrectPhi
    (
        U,
        phi,
        p_rgh,
        dimensionedScalar("rAUf", dimTime/rho.dimensions(), 1),
        closureSrc,   // zero at a fresh start
        pimple
    );

    #include "continuityErrs.H"
}
