

if (laser.powderSim())
{
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


Marangoni = condensate*epsilon1*mixture.MarangoniForce(T);

if (mthd.valid())
{
    mthd->updateLorentzForce();
}

const volVectorField& JcrossB =
    mthd.valid()
  ? mthd->JcrossB()
  : zeroJcrossB;

// inertia and gravity use the ledger density rhoM, the mass rhoPhi transports
fvVectorMatrix UEqn
(
    fvm::ddt(rhoM, U)
  + fvm::div(mixture.rhoPhi(), U)
  + turbulence->divDevRhoReff(U)
  + fvm::Sp(DC, U)
  + Marangoni
  - JcrossB
);

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
                mixture.surfaceTensionForce()
              - ghf*fvc::snGrad(rhoM)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    K = 0.5*magSqr(U);
}
