{
    TSolidus = alpha1*Tsolidus1 + (1.0 - alpha1)*Tsolidus2;
    TLiquidus = alpha1*Tliquidus1 + (1.0 - alpha1)*Tliquidus2;
    beta = alpha1*beta1 + (1.0 - alpha1)*beta2;
    rhok = (1.0 - max(epsilon1*(beta)*(T - TSolidus), 0.0));

    scalarField& alpha1I = alpha1;
    scalarField& TSolidusI = TSolidus;
    scalarField& TLiquidusI = TLiquidus;
    const scalar Tsolidus1Value = Tsolidus1.value();
    const scalar Tliquidus1Value = Tliquidus1.value();
    forAll(alpha1I, celli)
    {
        if (alpha1I[celli] > 0.05)
        {
            TSolidusI[celli] = Tsolidus1Value;
            TLiquidusI[celli] = Tliquidus1Value;
        }
    }

    TSolidus.correctBoundaryConditions();
    TLiquidus.correctBoundaryConditions();
    beta.correctBoundaryConditions();  // ADDED
    rhok.correctBoundaryConditions();  // ADDED

    #include "updateKappaCp.H"

    // Update the cell size fields
    const faceList & ff = mesh.faces();
    const pointField & pp = mesh.points();
    const vectorField& CI = mesh.C();
    scalarField& yDimI = laser.yDim();
    forAll(CI, celli)
    {
        const cell& cc = mesh.cells()[celli];
        labelList pLabels(cc.labels(ff));
        pointField pLocal(pLabels.size(), vector::zero);

        forAll(pLabels, pointi)
        {
            pLocal[pointi] = pp[pLabels[pointi]];
        }

        yDimI[celli] =
            Foam::max(pLocal & vector(0, 1, 0))
          - Foam::min(pLocal & vector(0, 1, 0));
    }

    laser.yDim().correctBoundaryConditions();

    metalaverage = fvc::average(alpha1);
    metalaverage.correctBoundaryConditions();  // ADDED

    alpha_filtered = volScalarField("alpha_filtered", alpha1);
    alpha_filtered.correctBoundaryConditions();
    scalarField& alpha_filteredI = alpha_filtered;
    // alpha_filteredI.correctBoundaryConditions();
    forAll(CI, celli)
    {
        if (alpha1I[celli] > 0.99)
        {
            alpha_filteredI[celli] = 1.0;
        }
        else
        {
            if (alpha1I[celli] < 0.01)
            {
                alpha_filteredI[celli] = 0.0;
            }
            else
            {
                alpha_filteredI[celli] = alpha1I[celli];
            }
        }
    }

    // Access the boundary fields
    volScalarField::Boundary& alpha_filteredBf = alpha_filtered.boundaryFieldRef();
    const volScalarField::Boundary& alpha1Bf = alpha1.boundaryField();

    // Loop through all boundary patches
    forAll(alpha_filteredBf, patchi)
    {
        // Get references to the patch fields
        scalarField& alpha_filteredPatch = alpha_filteredBf[patchi];
        const scalarField& alpha1Patch = alpha1Bf[patchi];

        // Apply the same logic to each face in the patch
        forAll(alpha_filteredPatch, facei)
        {
            if (alpha1Patch[facei] > 0.99)
            {
                alpha_filteredPatch[facei] = 1.0;
            }
            else if (alpha1Patch[facei] < 0.01)
            {
                alpha_filteredPatch[facei] = 0.0;
            }
            else
            {
                alpha_filteredPatch[facei] = alpha1Patch[facei];
            }
        }
    }

    alpha_filtered.correctBoundaryConditions();

    alpha_smoothed = fvc::average(alpha1);

    LatentHeat =
        alpha_filtered*LatentHeat1 + (1.0 - alpha_filtered)*LatentHeat2;

    // gradepsilon1 = fvc::grad(epsilon1);
    // gradepsilon1.correctBoundaryConditions();

    {
        const volVectorField gradAlphaFiltered(fvc::grad(alpha_filtered));
        n_filtered = gradAlphaFiltered/(mag(gradAlphaFiltered) + deltaN);
    }

    e1temp = fvc::average(epsilon1);

    const scalarField& e1tempI = e1temp;
    const scalarField& epsilon1I = epsilon1;
    scalarField& epsilon1maskI = epsilon1mask;
    forAll(CI, celli)
    {
        if (e1tempI[celli] <= 0.95)
        {
            epsilon1maskI[celli] = 0.0;
        }
        else
        {
            epsilon1maskI[celli] = epsilon1I[celli];
        }
    }

    epsilon1mask.correctBoundaryConditions();

    // e2temp = fvc::average(epsilon1mask);
    // epsilon1mask2 *= 0.0;

    // const scalarField& e2tempI = e2temp;
    // scalarField& epsilon1mask2I = epsilon1mask2;
    // forAll(CI, celli)
    // {
    //     if (e2tempI[celli] <= 0.95)
    //     {
    //         epsilon1mask2I[celli] = 0.0;
    //     }
    //     else
    //     {
    //         epsilon1mask2I[celli] = epsilon1maskI[celli];
    //     }

    // }

    // epsilon1mask2.correctBoundaryConditions();

    if (laser.powderSim())
    {
        gh = epsilon1mask*(g & mesh.C());
        ghf = fvc::interpolate(epsilon1mask)*(g & mesh.Cf());
    }
    else{
        gh = (g & mesh.C());
        ghf = (g & mesh.Cf());
    }
}
