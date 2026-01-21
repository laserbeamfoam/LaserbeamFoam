MRF.correctBoundaryVelocity(U);

gradAlpha = fvc::grad(alpha1);
gradAlpha.correctBoundaryConditions();  // ADDED - CRITICAL!

nHatM = gradAlpha/(mag(gradAlpha) + deltaN);
nHatM.correctBoundaryConditions();  // ADDED

Marangoni =
    epsilon1*Marangoni_Constant*(gradT - nHatM*(nHatM & gradT))*mag(gradAlpha);

const scalar Tliquidus1Value = Tliquidus1.value();
const scalarField& TI = T;
vectorField& MarangoniI = Marangoni;
forAll(MarangoniI, celli)
{
    if(TI[celli] < Tliquidus1Value)
    {
        MarangoniI[celli] = vector::zero;
    }
}

Marangoni.correctBoundaryConditions();

if (damperSwitch)
{
    // To account for high density difference
    damper = 2.0*rho/(rho1 + rho2);
    damper.correctBoundaryConditions();  // ADDED
}

pVap = 0.54*p0*Foam::exp(LatentHeatVap*Mm*((T - Tvap)/(R*T*Tvap)));
pVap.correctBoundaryConditions();  // ADDED


    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
            + fvm::Sp(solidDampingCoeff, U)
      - Marangoni*damper
     ==
                fvOptions(rho, U)
            + solidDampingCoeff*fvc::ddt(D)
            + fvc::div(sigmaElast)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()*fvc::interpolate(damper)
                    + (fvc::interpolate(pVap)
                        *fvc::snGrad(alpha1)
                        *fvc::interpolate(damper))
                  - ghf*fvc::snGrad(rho*rhok)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
