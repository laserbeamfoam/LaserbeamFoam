{
    alphasum *= 0.0;
    condensate = 1.0;

    TSolidus *= 0.0;
    TLiquidus *= 0.0;
    LatentHeat *= 0.0;
    rhoLatentHeat *= 0.0;   // sum_cond c_i*L_i [J/m^3], see createFields.H

    const dimensionedScalar cZeroLH("cZeroLH", dimDensity, Zero);
    label phaseiLH = 0;     // index into mixture.c(), phases() iteration order

    // electrical resistivity of condensed mixture/
    electrical_resistivity *= 0.0;

    if (mthd.valid())
    {
        mthd->resetTransportProperties();
    }



     volScalarField vapourSum
    (
        IOobject
        (
            "vapourSum",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("vapourSum", dimless, 0.0)
    );


    forAllConstIter(PtrDictionary<phaseModel>, mixture.phases(), iter)
    {
        const phaseModel& phase1 = iter();
        const volScalarField& alpha1 = phase1;

        // the phase's thermo is this dictionary (thermophysicalProperties.<phase>), already in memory
        const dictionary& phasedictionary = phase1.thermo();

        dimensionedScalar phaseTSolidus
        (
            "phaseTSolidus",
            dimensionSet(0, 0, 0, 1, 0),
            readScalar(phasedictionary.lookup("TSolidus"))
        );

        dimensionedScalar phaseTLiquidus
        (
            "phaseTLiquidus",
            dimensionSet(0, 0, 0, 1, 0),
            readScalar(phasedictionary.lookup("TLiquidus"))
        );

        dimensionedScalar phaseLH
        (
            "phaseLH",
            dimensionSet(0, 2, -2, 0, 0),
            readScalar(phasedictionary.lookup("LatentHeat"))
        );

        dimensionedScalar phaseconductivity
        (
            "phaseconductivity",
            dimensionSet(-1, -3, 3, 0, 0),
            readScalar(phasedictionary.lookup("elec_conductivity"))
        );

        dimensionedScalar phasepermeability
        (
            "phasepermeability",
            dimensionSet(1, 1, -2, 0, 0),
            readScalar(phasedictionary.lookup("magnetic_permeability"))
        );


        // Melt properties are a condensed-phase quantity.
        if (phase1.isCondensed())
        {
            TSolidus   += alpha1*phaseTSolidus;
            TLiquidus  += alpha1*phaseTLiquidus;
            LatentHeat += alpha1*phaseLH;

            // L is per unit mass, so the volume-additive latent heat is c_i*L with the ledger mass: the
            // field the energy equation uses (LatentHeat above is the volume average, output only)
            volScalarField cLH(max(mixture.c()[phaseiLH], cZeroLH));
            if (ledgerFloorEOS)
            {
                cLH = max(cLH, alpha1*phase1.thermo().rho());   // same guard as rhoCv
            }
            rhoLatentHeat += cLH*phaseLH;
        }

        if (mthd.valid())
        {
            mthd->accumulateTransportProperties
            (
                alpha1,
                phaseconductivity,
                phasepermeability
            );
        }

        if (phase1.isGaseous())
        {
            condensate -= alpha1;

                       if (phase1.name().find("vapour") != string::npos)
            {
                vapourSum += alpha1;
            }



        }

        if (phase1.isCondensed())
        {
            electrical_resistivity += alpha1/phaseconductivity;
        }




        alphasum += alpha1;
        ++phaseiLH;
    }

    condensate = max(min(condensate, 1.0), 0.0);

    electrical_resistivity /= max(condensate, SMALL);


    electrical_resistivity.correctBoundaryConditions();

    if (mthd.valid())
    {
        mthd->finaliseTransportProperties();
    }

    LatentHeat.correctBoundaryConditions();
    rhoLatentHeat.correctBoundaryConditions();
    TLiquidus.correctBoundaryConditions();
    TSolidus.correctBoundaryConditions();


    // Cell extent in y for the laser, rebuilt only when the mesh has changed (cellSizeStale: set in
    // createFields.H and after a mesh update)
    if (cellSizeStale)
    {
        const faceList& ff = mesh.faces();
        const pointField& pp = mesh.points();
        const cellList& cellFaces = mesh.cells();
        scalarField& yDimI = laser.yDim();

        forAll(cellFaces, celli)
        {
            scalar yMin = GREAT;
            scalar yMax = -GREAT;

            for (const label facei : cellFaces[celli])
            {
                for (const label pointi : ff[facei])
                {
                    const scalar y = pp[pointi].y();
                    yMin = Foam::min(yMin, y);
                    yMax = Foam::max(yMax, y);
                }
            }

            yDimI[celli] = yMax - yMin;
        }

        laser.yDim().correctBoundaryConditions();

        cellSizeStale = false;
    }

    metalaverage = fvc::average(condensate);
    condensateFiltered = max(min(condensate, scalar(1)), scalar(0));
    condensateFiltered.correctBoundaryConditions();



        {
        const dimensionedScalar betaBand("betaBand", dimless, 0.001);
        const dimensionedScalar vapRefineThresh
        (
            "vapRefineThresh", dimless, 0.001
        );

        const volScalarField vapourAverage(fvc::average(vapourSum));


        const dimensionedScalar coarsenGuard("coarsenGuard", dimless, 0.02);
        const volScalarField minorityFrac
        (
            min(condensateFiltered, scalar(1) - condensateFiltered)
        );

        refineIndicator =
            max
            (
                max
                (
                    min
                    (
                        min(metalaverage, scalar(1) - metalaverage)/betaBand,
                        scalar(1)
                    ),
                    min(vapourAverage/vapRefineThresh, scalar(1))
                ),
                min(minorityFrac/coarsenGuard, scalar(0.45))
            );
        refineIndicator.correctBoundaryConditions();
    }





    const volVectorField gradcondensate(fvc::grad(condensateFiltered));

    n_filtered = gradcondensate/(mag(gradcondensate) + deltaN);
    n_filtered.correctBoundaryConditions();
}
