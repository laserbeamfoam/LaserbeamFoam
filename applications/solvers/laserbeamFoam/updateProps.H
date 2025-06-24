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
    laser.refineFlag() *= 0.0;
    laser.refineFlag() += metalaverage;
    laser.refineFlag().correctBoundaryConditions();

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

    alpha_filtered.correctBoundaryConditions();

    LatentHeat =
        alpha_filtered*LatentHeat1 + (1.0 - alpha_filtered)*LatentHeat2;

    // gradepsilon1 = fvc::grad(epsilon1);
    // gradepsilon1.correctBoundaryConditions();

    n_filtered = fvc::grad(alpha_filtered)/((mag(fvc::grad(alpha_filtered)) + deltaN));

    e1temp = fvc::average(epsilon1);
    e1temp.correctBoundaryConditions();

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
