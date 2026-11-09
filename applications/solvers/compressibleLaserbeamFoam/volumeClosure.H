// Volume closure after TEqn, the pEqn's only volume source: s = sum_i c_i/rho_i(p, T) must stay 1 (the
// phase-norming of Zenz et al. 2024, Eqs. 53-62, https://doi.org/10.1016/j.compfluid.2023.106109).
// chi*(s - 1)/dt with chi = rhoCv/(rhoCv + p*ds/dT), the adiabatic share: 1/gamma in a gas, 1 in a liquid
volScalarField chiClosure
(
    IOobject
    (
        "chiClosure",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, 1.0)
);

closureSrc = mixture.volumeClosure(chiClosure);
closureSrc.correctBoundaryConditions();

// the p_rgh s was measured at, with this step's rhoM, so p_rghS + rhoM*gh == p exactly
const volScalarField p_rghS("p_rghS", p - rhoM*gh);
