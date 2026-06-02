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



TRHS=          LatentHeat*((
                fvc::ddt(rho,epsilon1)
              + fvc::div(mixture.rhoPhi(), epsilon1)
              )


);

TRHS.correctBoundaryConditions();


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
