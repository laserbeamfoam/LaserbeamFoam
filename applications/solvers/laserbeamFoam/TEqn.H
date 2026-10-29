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

    // apparentCapacity: in the cells of the mushy range the liquid fraction
    // is a function of T, epsilon1 = (T - TSolidus)/(TLiquidus - TSolidus), so
    // its latent heat is linearised about the point of that line for the
    // current liquid fraction, Tcorr = TSolidus + (TLiquidus - TSolidus)
    // *epsilon1, and taken implicitly, rho*L/(TLiquidus - TSolidus)
    // *(T - Tcorr)/deltaT (a Newton linearisation), and the correction loop
    // converges in a few corrections. The converged solution is the same as
    // with sourceBased
    const bool newtonLatentHeat = (latentHeatMethod == "apparentCapacity");

    // The coefficient of the linearised latent heat, only with
    // apparentCapacity
    autoPtr<volScalarField> latentCoeffPtr;

    if (newtonLatentHeat)
    {
        latentCoeffPtr.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "latentCoeff",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    rho.dimensions()*LatentHeat.dimensions()
                   /dimTemperature/dimTime,
                    Zero
                )
            )
        );
    }

    // Whether the latent heat of a cell is linearised, from its liquid
    // fraction and temperature: a cell in the mushy range, or a solid or
    // liquid cell whose temperature has entered the range
    const auto inMushyRange = [&](const label celli)
    {
        const scalar eps = epsilon1[celli];

        return
            TLiquidus[celli] - TSolidus[celli] > SMALL
         && LatentHeat[celli] > SMALL
         && (
                (eps > 0 && eps < 1)
             || (eps <= 0 && T[celli] > TSolidus[celli])
             || (eps >= 1 && T[celli] < TLiquidus[celli])
            );
    };

    do
    {
        epsilon1.storePrevIter();

        if (newtonLatentHeat)
        {
            // The cells in the mushy range, from the liquid fraction and the
            // temperature of the previous correction (a cell that leaves the
            // range is treated as solid or liquid in the next correction)
            scalarField& latentCoeffI = latentCoeffPtr->primitiveFieldRef();
            const scalarField& TSolidusI = TSolidus;
            const scalarField& TLiquidusI = TLiquidus;
            const scalarField& rhoI = rho;
            const scalarField& LatentHeatI = LatentHeat;
            const scalar rDeltaT = 1.0/runTime.deltaTValue();

            forAll(latentCoeffI, celli)
            {
                const scalar range = TLiquidusI[celli] - TSolidusI[celli];

                latentCoeffI[celli] =
                    inMushyRange(celli)
                  ? rhoI[celli]*LatentHeatI[celli]/range*rDeltaT
                  : 0;
            }

            Tcorr = (TLiquidus - TSolidus)*epsilon1 + TSolidus;
            Tcorr.correctBoundaryConditions();
        }
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
    const volScalarField& Joulesource =
        mthd.valid()
      ? mthd->Joulesource()
      : zeroJoulesource;

        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp,T)
          + fvm::div(rhophicp, T)
          - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhophicp), T)
          - fvm::laplacian(kappa, T)
          - fvm::Sp(Qrad_inpl*gradAlphaDamper, T)
          ==
            laser.deposition()
          + Joulesource
          + (Qrad_expl - Qv)*gradAlphaDamper
          - TRHS
        );

        if (newtonLatentHeat)
        {
            // rho*L/(TLiquidus - TSolidus)*(T - Tcorr)/deltaT: implicit in T
            const volScalarField& latentCoeff = latentCoeffPtr();
            TEqn += fvm::Sp(latentCoeff, T);
            TEqn -= latentCoeff*Tcorr;
        }

        TEqn.relax();

        TEqn.solve();

        T.correctBoundaryConditions();

        // T.relax();

        if (newtonLatentHeat)
        {
            // Liquid fraction from the new temperature: the Newton update in
            // the cells of the mushy range. A solid or liquid cell keeps its
            // liquid fraction; if its new temperature has entered the mushy
            // range, it is classified as mushy in the next correction and
            // linearised from the edge of the range (Tcorr = TSolidus or
            // TLiquidus), rather than set to the equilibrium value of a
            // temperature that may have overshot the whole range
            scalarField& epsilon1I = epsilon1.primitiveFieldRef();
            const scalarField& latentCoeffI = latentCoeffPtr();
            const scalarField& TcorrI = Tcorr;
            const scalarField& TI = T;
            const scalarField& TSolidusI = TSolidus;
            const scalarField& TLiquidusI = TLiquidus;

            forAll(epsilon1I, celli)
            {
                const scalar range = TLiquidusI[celli] - TSolidusI[celli];

                if (range <= SMALL)
                {
                    continue;
                }

                if (latentCoeffI[celli] > 0)
                {
                    const scalar eps =
                        epsilon1I[celli] + (TI[celli] - TcorrI[celli])/range;

                    epsilon1I[celli] = max(min(eps, scalar(1)), scalar(0));
                }
            }
        }
        else
        {
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
        }

        epsilon1.correctBoundaryConditions();

        residual =
            mag
            (
                epsilon1.primitiveField()
              - epsilon1.prevIter().primitiveField()
            );

        if (newtonLatentHeat)
        {
            // A solid or liquid cell whose new temperature has entered the
            // mushy range keeps its liquid fraction until the next correction:
            // its residual is how far its temperature has entered the range,
            // so that the loop does not stop before its latent heat is taken
            const scalarField& latentCoeffI = latentCoeffPtr();
            const scalarField& epsilon1I = epsilon1;
            const scalarField& TI = T;
            const scalarField& TSolidusI = TSolidus;
            const scalarField& TLiquidusI = TLiquidus;

            forAll(residual, celli)
            {
                if (latentCoeffI[celli] <= 0 && inMushyRange(celli))
                {
                    const scalar range = TLiquidusI[celli] - TSolidusI[celli];

                    residual[celli] =
                        min
                        (
                            epsilon1I[celli] <= 0
                          ? (TI[celli] - TSolidusI[celli])/range
                          : (TLiquidusI[celli] - TI[celli])/range,
                            scalar(1)
                        );
                }
            }
        }

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
