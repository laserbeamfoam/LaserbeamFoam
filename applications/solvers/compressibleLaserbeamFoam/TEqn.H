{

  //https://www.cfd-online.com/Forums/openfoam-solving/251323-teqn-compressibleinterfoam.html





    int iter = 0;


    dimensionedScalar residual("residual",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));
    dimensionedScalar meanResidual("meanResidual",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));
    dimensionedScalar residual2("residual2",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));
    dimensionedScalar meanResidual2("meanResidual2",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));





    TRHS*=0.0;//=TRHSzero;

    do{
iter++;
epsilon1.storePrevIter();

TRHS=         LatentHeat*(
              (
                fvc::ddt(rho,epsilon1)
              + fvc::div(mixture.rhoPhi(), epsilon1)
              )
);

TRHS.correctBoundaryConditions();

// RB: leaving this temporarily in case anyone is checking solver-restarting
// during review.
if (runTime.writeTime() && debugFlag)
{
    volScalarField& rhoCopyPreUpdate = rhoCopyPreUpdatePtr();

    volScalarField ddt_term_error
    (
        IOobject
        (
            "ERROR_ddt_rho_euler_scheme",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
            fvc::ddt(rho, epsilon1)
        - ( epsilon1*(rho - rho.oldTime())/mesh.time().deltaT())
    );

    volScalarField ddt_term
    (
        IOobject
        (
            "TRHS_ddt_rho_epsilon1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::ddt(rho,epsilon1)
    );

    volScalarField delta_rho_rhoOldTime
    (
        IOobject
        (
            "Diff_rho_rhoOldTime",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
           (rho - rho.oldTime())
    );

    volScalarField delta_rho_rhoCopy
    (
        IOobject
        (
            "Diff_rho_rhoCopyPreUpdate",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
           (rho - rhoCopyPreUpdate)
    );

    volScalarField delta_rhocopy_in_time
    (
        IOobject
        (
            "Diff_rhoCopyPreUpdate_rhoCopyPreUpdate.oldTime",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
           (rhoCopyPreUpdate - rhoCopyPreUpdate.oldTime())
    );

    volScalarField ddt_term_rho_only
    (
        IOobject
        (
            "TRHS_ddt_rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::ddt(rho)
    );

    volScalarField div_term
    (
        IOobject
        (
            "TRHS_div_rhoPhi_epsilon1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::div(mixture.rhoPhi(), epsilon1)
    );

    const volScalarField rhoOld_check
    (
        IOobject
        (
            "rhoOldTime_during_run",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho.oldTime()
    );

    delta_rhocopy_in_time.write();
    rhoOld_check.write();

    ddt_term.write();
    ddt_term_rho_only.write();
    div_term.write();
    ddt_term_error.write();

    delta_rho_rhoOldTime.write();
    rhoCopyPreUpdate.write();
    delta_rho_rhoCopy.write();
}

    const volScalarField& Joulesource =
        mthd.valid()
      ? mthd->Joulesource()
      : zeroJoulesource;

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T)
      + fvm::div(mixture.rhoPhi(), T)
      - fvm::laplacian(mixture.alphaEff(turbulence->mut()), T)
      + (
            fvc::div(fvc::absolute(phi, U), p)
          + fvc::ddt(rho, K) + fvc::div(mixture.rhoPhi(), K)
        )*mixture.rCv()

        ==
          mixture.rCv()*laser.deposition()
        + mixture.rCv()*Joulesource
        - mixture.rCv()*TRHS
        -mass_dot
    );

    TEqn.relax();
    TEqn.solve();


Tcorr=(TLiquidus-TSolidus)*epsilon1+TSolidus;
Tcorr.correctBoundaryConditions();


epsilon1=max(min(epsilon1+(epsilonRel*mixture.Cv()/LatentHeat)*(T-Tcorr),scalar(1)),scalar(0));
epsilon1.correctBoundaryConditions();


residual = gMax(mag(epsilon1.primitiveField()-epsilon1.prevIter().primitiveField()));
meanResidual = gSum(mag(epsilon1.primitiveField()-epsilon1.prevIter().primitiveField())*mesh.V())/gSum(mesh.V().field());


        Info << "Correcting epsilon1, mean residual = " << meanResidual.value()
             << ", max residual = " << residual.value()
             << endl;


}while((iter < minTCorr || (residual > epsilonTol)) && iter <= maxTCorr );



    mixture.correct();
}
