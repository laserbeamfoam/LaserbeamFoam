{


    int iter = 0;
    scalarField residual(epsilon1.primitiveField().size(), 1.0);
    scalar maxResidual = 1.0;
    scalar meanResidual = 1.0;

    TRHS *= 0.0;


    #include "updateRhoCv.H"


    rhoCv.oldTime() = rho.oldTime()*CvMass.oldTime();

    // Cv-weighted mass flux (laserbeamFoam: rhophicp = interpolate(cp)*rhoPhi)
    const surfaceScalarField rhoCvPhi
    (
        "rhoCvPhi",
        fvc::interpolate(CvMass)*mixture.rhoPhi()
    );


    const volScalarField kappaEff
    (
        "kappaEff",
        mixture.kappaEff(turbulence->mut())
    );


    const volScalarField contErrCv
    (
        "contErrCv",
        fvc::ddt(rhoCv) + fvc::div(rhoCvPhi)
    );
    const volScalarField contErrRho
    (
        "contErrRho",
        fvc::ddt(rho) + fvc::div(mixture.rhoPhi())
    );

    // ---- melting / solidification -----------------------------------------

    const dimensionedScalar cMinMelt   ("cMinMelt",    dimless, 1e-4); // div guard
    const dimensionedScalar cThreshMelt("cThreshMelt", dimless, 1e-2); // mask ramp

    const volScalarField rCond(1.0/max(condensate, cMinMelt));
    const volScalarField TSolidusCond  (TSolidus*rCond);    // -> metal T_sol
    const volScalarField TLiquidusCond (TLiquidus*rCond);   // -> metal T_liq


    const volScalarField LatentHeatMass
    (
        "LatentHeatMass",
        rhoLatentHeat/max(rho, dimensionedScalar("rhoSmallT", dimDensity, SMALL))
    );

    const dimensionedScalar rhoLfloor
    (
        "rhoLfloor", rhoLatentHeat.dimensions(), 1.0
    );


    const volScalarField meltCoeff
    (
        epsilonRel*rhoCvCond/max(rhoLatentHeat, rhoLfloor)
    );

    // 1 where metal present, ramps to 0 in the gas
    const volScalarField meltMask
    (
        min(max(condensate/cThreshMelt, scalar(0)), scalar(1))
    );


    boolList e1MaskCell(mesh.nCells(), false);
    {
        const wordList e1MaskNames
        (
            mesh.lookupObject<IOdictionary>("thermophysicalProperties")
                .lookupOrDefault<wordList>("phaseChangeMaskPatches", wordList())
        );

        for (const word& pn : e1MaskNames)
        {
            const label patchi = mesh.boundaryMesh().findPatchID(pn);
            if (patchi < 0) continue;   // typo already warned about in solveAlphas

            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            forAll(fc, bfacei)
            {
                e1MaskCell[fc[bfacei]] = true;
            }
        }
    }


    const bool newtonLatentHeat = (latentHeatMethod == "apparentCapacity");
    const scalar rDeltaT = 1.0/runTime.deltaTValue();

    volScalarField latentCoeff
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
            rhoLatentHeat.dimensions()/dimTemperature/dimTime,
            Zero
        )
    );

    const scalarField& TSolidusCondI  = TSolidusCond.primitiveField();
    const scalarField& TLiquidusCondI = TLiquidusCond.primitiveField();
    const scalarField& rhoLatentHeatI = rhoLatentHeat.primitiveField();
    const scalarField& meltMaskI      = meltMask.primitiveField();
    const scalarField& meltCoeffI     = meltCoeff.primitiveField();


    const auto inMushyRange = [&](const label celli)
    {
        const scalar eps = epsilon1[celli];
        const scalar Tc  = T[celli];

        return
            meltMaskI[celli] >= 1.0 - SMALL
         && !e1MaskCell[celli]
         && TLiquidusCondI[celli] - TSolidusCondI[celli] > SMALL
         && rhoLatentHeatI[celli] > SMALL
         && (
                (eps > 0 && eps < 1)
             || (eps <= 0 && Tc > TSolidusCondI[celli])
             || (eps >= 1 && Tc < TLiquidusCondI[celli])
            );
    };

    label nNewton = 0;

    do
    {
        iter++;
        epsilon1.storePrevIter();

        // linearisation point of this correction (both methods)
        Tcorr = (TLiquidusCond - TSolidusCond)*epsilon1 + TSolidusCond;
        Tcorr.correctBoundaryConditions();

        if (newtonLatentHeat)
        {

            scalarField& latentCoeffI = latentCoeff.primitiveFieldRef();
            nNewton = 0;

            forAll(latentCoeffI, celli)
            {
                if (inMushyRange(celli))
                {
                    const scalar range =
                        TLiquidusCondI[celli] - TSolidusCondI[celli];

                    latentCoeffI[celli] = rhoLatentHeatI[celli]/range*rDeltaT;
                    ++nNewton;
                }
                else
                {
                    latentCoeffI[celli] = 0;
                }
            }
            reduce(nNewton, sumOp<label>());
        }

        TRHS =
            meltMask*LatentHeatMass
           *(
                fvc::ddt(rho, epsilon1)
              + fvc::div(mixture.rhoPhi(), epsilon1)
              - epsilon1*contErrRho
            );

        TRHS.correctBoundaryConditions();

        const volScalarField& Joulesource =
            mthd.valid()
          ? mthd->Joulesource()
          : zeroJoulesource;

        // every term in W/m^3
        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCv, T)
          + fvm::div(rhoCvPhi, T)
          - fvm::Sp(contErrCv, T)
          - fvm::laplacian(kappaEff, T)
          + (
                fvc::div(fvc::absolute(phi, U), p)
              + fvc::ddt(rho, K) + fvc::div(mixture.rhoPhi(), K)
            )
         ==
            laser.deposition()
          + Joulesource
          - TRHS
          - mass_dot
        );

        if (newtonLatentHeat)
        {
            // rhoLatentHeat/(Tl - Ts)*(T - Tcorr)/deltaT, implicit in T
            TEqn += fvm::Sp(latentCoeff, T);
            TEqn -= latentCoeff*Tcorr;
        }

        TEqn.relax();
        TEqn.solve();


        {
            const dimensionedScalar Tfloor("Tfloor", dimTemperature, 250.0);

            const scalar Tceil = pimple.dict().getOrDefault<scalar>("TMax", 1e6);
            scalarField& TI = T.primitiveFieldRef();
            label nClip = 0, nClipHi = 0;
            forAll(TI, celli)
            {
                if (TI[celli] < Tfloor.value()) { TI[celli] = Tfloor.value(); ++nClip; }
                if (TI[celli] > Tceil) { TI[celli] = Tceil; ++nClipHi; }
            }
            reduce(nClip, sumOp<label>());
            reduce(nClipHi, sumOp<label>());
            if (nClip) Info<< "TEqn: floored " << nClip << " cell(s)" << endl;
            if (nClipHi) Info<< "TEqn: ceilinged " << nClipHi << " cell(s)" << endl;
            T.correctBoundaryConditions();
        }


 
        scalarField newtonStep(epsilon1.primitiveField().size(), 0.0);
        {
            scalarField& e1I = epsilon1.primitiveFieldRef();
            const scalarField& TI = T.primitiveField();
            const scalarField& TcorrI = Tcorr.primitiveField();
            const scalarField& latentCoeffI = latentCoeff.primitiveField();

            forAll(e1I, celli)
            {
                const scalar epsOld = e1I[celli];
                scalar epsNew = epsOld;

                if (latentCoeffI[celli] > 0)
                {
                    const scalar range =
                        TLiquidusCondI[celli] - TSolidusCondI[celli];

                    const scalar epsUnclipped =
                        epsOld + (TI[celli] - TcorrI[celli])/range;

                    epsNew = Foam::max(Foam::min(epsUnclipped, scalar(1)), scalar(0));
                    newtonStep[celli] = Foam::mag(epsUnclipped - epsOld);
                }
                else
                {
                    epsNew = Foam::max
                    (
                        Foam::min
                        (
                            epsOld + meltCoeffI[celli]*(TI[celli] - TcorrI[celli]),
                            scalar(1)
                        ),
                        scalar(0)
                    );
                }

                e1I[celli] =
                    meltMaskI[celli]*epsNew + (scalar(1) - meltMaskI[celli]);
            }
        }


        {
            label nE1Mask = 0;
            scalarField& e1I = epsilon1.primitiveFieldRef();

            forAll(e1MaskCell, celli)
            {
                if (e1MaskCell[celli] && e1I[celli] < 1.0)
                {
                    e1I[celli] = 1.0;
                    ++nE1Mask;
                }
            }
            reduce(nE1Mask, sumOp<label>());
            if (nE1Mask)
            {
                Info<< "solidification masked in " << nE1Mask
                    << " boundary cell(s)" << endl;
            }
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
            const scalarField& latentCoeffI = latentCoeff.primitiveField();
            const scalarField& e1I = epsilon1.primitiveField();
            const scalarField& TI = T.primitiveField();

            forAll(residual, celli)
            {
                if (latentCoeffI[celli] > 0)
                {

                    residual[celli] = Foam::max(residual[celli], newtonStep[celli]);
                }
                else if (inMushyRange(celli))
                {
                    const scalar range =
                        TLiquidusCondI[celli] - TSolidusCondI[celli];

                    residual[celli] = Foam::max
                    (
                        residual[celli],
                        Foam::min
                        (
                            e1I[celli] <= 0
                          ? (TI[celli] - TSolidusCondI[celli])/range
                          : (TLiquidusCondI[celli] - TI[celli])/range,
                            scalar(1)
                        )
                    );
                }
            }
        }

        maxResidual = gMax(residual);
        meanResidual = gSum(residual*mesh.V())/gSum(mesh.V().field());

        Info<< "Correcting epsilon1 (" << latentHeatMethod;
        if (newtonLatentHeat)
        {
            Info<< ", " << nNewton << " cells linearised";
        }
        Info<< "), mean residual = " << meanResidual
            << ", max residual = " << maxResidual
            << endl;

    } while ((iter < minTCorr || (maxResidual > epsilonTol.value())) && iter <= maxTCorr);


    if (newtonLatentHeat)
    {
        const scalarField& latentCoeffI = latentCoeff.primitiveField();
        const scalarField& TI = T.primitiveField();
        const scalarField& TcorrI = Tcorr.primitiveField();
        const scalarField& e1I = epsilon1.primitiveField();
        const scalarField& e1PrevI = epsilon1.prevIter().primitiveField();
        const scalarField& V = mesh.V().field();

        scalar implicitW = 0.0;
        scalar realisedW = 0.0;
        forAll(latentCoeffI, celli)
        {
            if (latentCoeffI[celli] > 0)
            {
                implicitW += latentCoeffI[celli]*(TI[celli] - TcorrI[celli])*V[celli];
                realisedW +=
                    rhoLatentHeatI[celli]*(e1I[celli] - e1PrevI[celli])*rDeltaT*V[celli];
            }
        }
        reduce(implicitW, sumOp<scalar>());
        reduce(realisedW, sumOp<scalar>());

        Info<< "latent heat linearisation: implicit " << implicitW
            << " W, realised " << realisedW
            << " W, un-cancelled " << implicitW - realisedW << " W"
            << endl;
    }


    mixture.correct();
}