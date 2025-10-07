// zeroField Su;
zeroField Sp;
zeroField divU;


volScalarField Su
(
    IOobject
    (
        "Su",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alphaMeltSource  // All explicit alpha sources
);