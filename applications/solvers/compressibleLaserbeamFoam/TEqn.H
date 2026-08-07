{

  //https://www.cfd-online.com/Forums/openfoam-solving/251323-teqn-compressibleinterfoam.html





    int iter = 0;


    dimensionedScalar residual("residual",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));
    dimensionedScalar meanResidual("meanResidual",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));
    dimensionedScalar residual2("residual2",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));
    dimensionedScalar meanResidual2("meanResidual2",dimensionSet(0, 0, 0, 0, 0),scalar(1.0));





    TRHS*=0.0;//=TRHSzero;

    // no spurious mushy drag in gas
    const dimensionedScalar cMinMelt   ("cMinMelt",    dimless, 1e-4); // div guard
    const dimensionedScalar cThreshMelt("cThreshMelt", dimless, 1e-2); // mask ramp

    const volScalarField rCond(1.0/max(condensate, cMinMelt));
    const volScalarField LatentHeatCond(LatentHeat*rCond);  // -> metal L_f
    const volScalarField TSolidusCond  (TSolidus*rCond);    // -> metal T_sol
    const volScalarField TLiquidusCond (TLiquidus*rCond);   // -> metal T_liq

    // Apparent-heat-capacity coefficient
    // LatentHeat is now condensed-only (see update.H), so its  0 in
    // pure-gas cells where it used to be the gas dummy (~1).  
    const volScalarField meltCoeff
    (
        epsilonRel*mixture.Cv()
       /max(LatentHeatCond, dimensionedScalar("LHfloor", LatentHeatCond.dimensions(), 1.0))
    );

    // 1 where metal present, ramps to 0 in the gas
    const volScalarField meltMask
    (
        min(max(condensate/cThreshMelt, scalar(0)), scalar(1))
    );
    // ------------------------------------------------------------------------

    do{
iter++;
epsilon1.storePrevIter();



TRHS = meltMask*LatentHeat*
(
    fvc::ddt(rho, epsilon1)
  + fvc::div(mixture.rhoPhi(), epsilon1)
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
      - fvm::Sp(fvc::ddt(rho) + fvc::div(mixture.rhoPhi()), T)
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



 {
        const dimensionedScalar Tfloor("Tfloor", dimTemperature, 250.0);
        scalarField& TI = T.primitiveFieldRef();
        label nClip = 0;
        forAll(TI, celli)
            if (TI[celli] < Tfloor.value()) { TI[celli] = Tfloor.value(); ++nClip; }
        reduce(nClip, sumOp<label>());
        if (nClip) Info<< "TEqn: floored " << nClip << " cell(s)" << endl;
        T.correctBoundaryConditions();
    }


Tcorr = (TLiquidusCond - TSolidusCond)*epsilon1 + TSolidusCond;
Tcorr.correctBoundaryConditions();

// Masked, normalised melt-fraction update. In metal (meltMask=1) this is
// identical to the original; in (near-)gas (meltMask->0) epsilon1 -> 1.
epsilon1 =
    meltMask*max(min(epsilon1 + meltCoeff*(T - Tcorr), scalar(1)), scalar(0))
  + (scalar(1) - meltMask);
epsilon1.correctBoundaryConditions();


residual = gMax(mag(epsilon1.primitiveField()-epsilon1.prevIter().primitiveField()));
meanResidual = gSum(mag(epsilon1.primitiveField()-epsilon1.prevIter().primitiveField())*mesh.V())/gSum(mesh.V().field());


        Info << "Correcting epsilon1, mean residual = " << meanResidual.value()
             << ", max residual = " << residual.value()
             << endl;


}while((iter < minTCorr || (residual > epsilonTol)) && iter <= maxTCorr );



    mixture.correct();
}