{
    TSolidus = alpha1*Tsolidus1 + (1.0 - alpha1)*Tsolidus2;
    TLiquidus = alpha1*Tliquidus1 + (1.0 - alpha1)*Tliquidus2;

    if (mthd.valid())
    {
        mthd->correctIncompressibleProperties
        (
            alpha1,
            transportPropertiesMetal,
            transportPropertiesGas
        );
    }

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
    beta.correctBoundaryConditions();
    rhok.correctBoundaryConditions();

    #include "updateKappaCp.H"

    const vectorField& CI = mesh.C();
    metalaverage = fvc::average(alpha1);
    metalaverage.correctBoundaryConditions();
    alpha_smoothed = metalaverage;
    alpha_smoothed.correctBoundaryConditions();

    scalarField& alpha_filteredI = alpha_filtered;
    alpha_filteredI = alpha1I;

    forAll(CI, celli)
    {
        const scalar alpha1Value = alpha1I[celli];

        if (alpha1Value > 0.99)
        {
            alpha_filteredI[celli] = 1.0;
        }
        else if (alpha1Value < 0.01)
        {
            alpha_filteredI[celli] = 0.0;
        }
        else
        {
            alpha_filteredI[celli] = alpha1Value;
        }
    }

    volScalarField::Boundary& alpha_filteredBf = alpha_filtered.boundaryFieldRef();
    const volScalarField::Boundary& alpha1Bf = alpha1.boundaryField();

    forAll(alpha_filteredBf, patchi)
    {
        scalarField& alpha_filteredPatch = alpha_filteredBf[patchi];
        const scalarField& alpha1Patch = alpha1Bf[patchi];

        forAll(alpha_filteredPatch, facei)
        {
            const scalar alpha1Value = alpha1Patch[facei];

            if (alpha1Value > 0.99)
            {
                alpha_filteredPatch[facei] = 1.0;
            }
            else if (alpha1Value < 0.01)
            {
                alpha_filteredPatch[facei] = 0.0;
            }
            else
            {
                alpha_filteredPatch[facei] = alpha1Value;
            }
        }
    }

    alpha_filtered.correctBoundaryConditions();

    LatentHeat =
        alpha_filtered*LatentHeat1 + (1.0 - alpha_filtered)*LatentHeat2;

    {
        const volVectorField gradAlphaFiltered(fvc::grad(alpha_filtered));
        n_filtered = gradAlphaFiltered/(mag(gradAlphaFiltered) + deltaN);
    }
    n_filtered.correctBoundaryConditions();

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

    if (laser.powderSim())
    {
        gh = epsilon1mask*(g & mesh.C());
        ghf = fvc::interpolate(epsilon1mask)*(g & mesh.Cf());
    }
    else
    {
        gh = (g & mesh.C());
        ghf = (g & mesh.Cf());
    }
}
