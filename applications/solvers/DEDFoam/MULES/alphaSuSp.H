// ============================================================================
// MULES/alphaSuSp.H
//
// Defines the Su (explicit) and Sp (implicit) source coefficient fields used
// by the MULES alpha equation. The standard interFoam version uses
// `zeroField Su; zeroField Sp;`. Here Su is promoted to a real volScalarField
// that contains alphaMeltSource (the volumetric rate at which particle mass
// turns into metal, [1/s]), so the alpha equation sees the melt addition.
//
// Su units must match the residual of the alpha equation, i.e. [1/s]. That
// matches alphaMeltSource which is also [1/s].
// ============================================================================

// zeroField Su;   // Replaced by the volScalarField below
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