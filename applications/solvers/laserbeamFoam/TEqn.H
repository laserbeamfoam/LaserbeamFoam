{
    int iter = 0;
    scalarField residual(epsilon1.primitiveField().size(), 1.0);
    scalar maxResidual = 1.0;
    scalar meanResidual = 1.0;
    //scalar Tresidual = 1.0;

    TRHS *= 0.0;

    rhoCp = rho*cp;
    rhoCp.correctBoundaryConditions();  
    
    rhophicp = fvc::interpolate(cp)*rhoPhi;

    em = (epsilon1*emL) + ((1.0 - epsilon1)*emS);

    // const volTensorField gradU(fvc::grad(U));
    // volTensorField tau
    // (
    //     ((mixture.mu()*gradU) + (mixture.mu()*gradU.T()))
    //   - ((2.0/3.0)*mixture.mu()*(fvc::div(phi))*I)
    // );

    do
    {
        epsilon1.storePrevIter();
        T.storePrevIter();
        TRHS =
            LatentHeat*((fvc::ddt(rho, epsilon1) + fvc::div(rhoPhi, epsilon1)));
        
        TRHS.correctBoundaryConditions();  

        Qv =
            0.82*LatentHeatVap*Mm*p0
           *Foam::exp(LatentHeatVap*Mm*((T - Tvap)/(R*T*Tvap)))
           /Foam::pow(2*M_PI*Mm*R*T,0.5);
        
        Qv.correctBoundaryConditions();  

        const volScalarField Qrad_expl (SB*em*(Foam::pow(TRef,4) + 3*Foam::pow(T,4)));
        const volScalarField Qrad_inpl (-(4*SB*em*Foam::pow(T,3)));

        if (damperSwitch)
        {
            scalarField& thermalDamperI = thermalDamper;
            const scalarField& rhoCpI = rhoCp;
            const scalarField& TI = T;
            forAll(thermalDamperI, celli)
            {
                thermalDamperI[celli] =
                    2.0*rhoCpI[celli]
                   /(
                       (rho1.value()*polycp_m.value(TI[celli]))
                     + (rho2.value()*polycp_g.value(TI[celli]))
                   );
            }

            thermalDamper.correctBoundaryConditions();
        }

        const volScalarField gradAlphaDamper (mag(gradAlpha)*thermalDamper);
        
        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp,T)
          + fvm::div(rhophicp, T)
          - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhophicp), T)
          - fvm::laplacian(kappa, T)
          - fvm::Sp(Qrad_inpl*gradAlphaDamper, T)
          //- (tau && gradU) // Viscous dissipation
         ==
            laser.deposition()
          - Qv*mag(gradAlpha)*thermalDamper
          + (Qrad_expl - Qv)*gradAlphaDamper
          - TRHS
        );

        TEqn.relax();

        TEqn.solve();
        
        T.correctBoundaryConditions();  

        // T.relax();

        Tcorr = (TLiquidus - TSolidus)*epsilon1 + TSolidus;
        Tcorr.correctBoundaryConditions();  

        epsilon1 =
            max
            (
                min
                (
                    epsilon1 + (epsilonRel*cp/LatentHeat)*(T - Tcorr),
                    scalar(1)
                ),
                scalar(0)
            );
        
        epsilon1.correctBoundaryConditions();  

        residual =
            mag
            (
                epsilon1.primitiveField()
              - epsilon1.prevIter().primitiveField()
            );

        maxResidual = gMax(residual);

        meanResidual = gSum(residual*mesh.V())/gSum(mesh.V().field());

        //ddte1 = fvc::ddt(epsilon1);

        Info<< "Correcting epsilon1, mean residual = " << meanResidual
            << ", max residual = " << maxResidual
            << endl;
    }
    while
    (
        (
            ++iter < minTCorr
         // || (meanResidual > epsilonTol.value())
         || (maxResidual > epsilonTol.value())
        )
     && iter <= maxTCorr
    );

    gradT = fvc::grad(T);
    gradT.correctBoundaryConditions();  // ADDED - CRITICAL!
}