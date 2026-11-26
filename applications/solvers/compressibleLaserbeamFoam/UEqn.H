// enthalpy-porosity (Carman-Kozeny) drag C*(1 - eps)^2/(eps^3 + b) holding the solid and the mush:
// Voller & Prakash 1987, https://doi.org/10.1016/0017-9310(87)90317-6 (C, b: fvSolution MELTING)
if (laser.powderSim())
{
    // powder: the drag is extended by one cell around the solid (as laserbeamFoam does), which damps
    // the spurious currents of the particles' staircased surfaces. Fluid cells with any solid
    // neighbour count as solid
    {
        const volScalarField epsilon1Avg(fvc::average(epsilon1));
        const scalarField& avgI = epsilon1Avg.primitiveField();
        const scalarField& e1I = epsilon1.primitiveField();
        scalarField& maskI = epsilon1mask.primitiveFieldRef();

        forAll(maskI, celli)
        {
            maskI[celli] = (avgI[celli] <= 0.95) ? 0.0 : e1I[celli];
        }
        epsilon1mask.correctBoundaryConditions();
    }

    DC =
        DarcyConstantlarge*Foam::pow((1.0-epsilon1mask),2)
       /(Foam::pow(epsilon1mask, 3) + DarcyConstantsmall);
}
else
{
    DC =
        DarcyConstantlarge*Foam::pow((1.0 - epsilon1), 2)
       /(Foam::pow(epsilon1, 3) + DarcyConstantsmall);
}

// thermocapillary stress dsigma/dT*grad_s(T)*|grad(beta)| put on the melt side by 2*beta, which keeps
// its integral across the surface at 1 (beta alone gave 1/2): the rho_l >> rho_g limit of the
// density-scaled delta 2*rho/(rho_l + rho_g), https://doi.org/10.1016/0021-9991(92)90240-Y
Marangoni = 2.0*condensate*epsilon1*mixture.MarangoniForce(T);

if (mthd.valid())
{
    mthd->updateLorentzForce();
}

const volVectorField& JcrossB =
    mthd.valid()
  ? mthd->JcrossB()
  : zeroJcrossB;

// the Stefan flow of this step's condensation (Upc, phiJump), once per time step
if (pimple.firstIter())
{
    #include "stefanFlow.H"
}

// inertia and gravity use the ledger density rhoM, the mass rhoPhi transports
fvVectorMatrix UEqn
(
    fvm::ddt(rhoM, U)
  + fvm::div(mixture.rhoPhi(), U)
  + turbulence->divDevRhoReff(U)
);

// viscous force of the Stefan flow: (4/3)*mu*grad(div(Upc)), zero off the surface, the spurious
// stress on it. Subtracted, so the viscous term acts on U - Upc (stefanFlow.H). Nothing to remove in
// a step without condensation (Upc is zero)
if (stefanFlowActive)
{
    // Upc under U's name (not registered), so its viscous operator is evaluated with U's schemes
    volVectorField UpcAsU
    (
        IOobject
        (
            U.name(),
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        Upc
    );

    const volVectorField viscousStefan
    (
        "viscousStefan",
        turbulence->divDevRhoReff(UpcAsU) & UpcAsU
    );

    UEqn -= viscousStefan;
}

UEqn += fvm::Sp(DC, U);
UEqn += Marangoni;
UEqn -= JcrossB;
UEqn -= mixture.phaseChangeMomentum();   // momentum the phase-change mass carries (solveAlphas)

if (ledgerFloorEOS)
{
    // guarded cells hold EOS mass rhoPhi never moved: rhoM*DU/Dt form, same as TEqn's Sp(contErrCv)
    const volScalarField contErrU(fvc::ddt(rhoM) + fvc::div(mixture.rhoPhi()));
    UEqn -= fvm::Sp(contErrU, U);
}


UEqn.relax();

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
                stf
              - ghf*fvc::snGrad(rhoM)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    K = 0.5*magSqr(U);
}
