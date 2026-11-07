// Volume closure on the post-TEqn state, every outer corrector: the pEqn's only volume source.
// s = sum_i c_i/rho_i(p, T): heating, phase change and transport all show up here once.
// chi*(s - 1)/dt with chi = rhoCv/(rhoCv + p*ds/dT): an outflow of chi*(s - 1) plus its own
// p*dV cooling removes the excess exactly (chi = 1/gamma in a gas, 1 in the liquid).
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
