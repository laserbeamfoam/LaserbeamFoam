{
    const scalarField& alpha1I = alpha1;
    const scalarField& massFractionI = massFraction;
    const scalarField& TI = T;
    scalarField& kappaI = kappa;
    scalarField& cpI = cp;
    scalarField& resistivityI = electrical_resistivity;
    forAll(kappaI, celli)
    {
        const scalar a = max(min(alpha1I[celli], 1.0), 0.0);
        const scalar Ti = TI[celli];

        kappaI[celli] =
            a*kappaTable_m.value(Ti)
          + (1.0 - a)*kappaTable_g.value(Ti);

        // cp is specific (per unit mass), so it is mixed with the mass
        // fraction so that rho*cp is volume-additive
        cpI[celli] =
            massFractionI[celli]*cpTable_m.value(Ti)
          + (1.0 - massFractionI[celli])*cpTable_g.value(Ti);

        // the laser evaluates Fresnel in the interface cells so they have to
        // look like metal: metal resistivity everywhere, not blended
        resistivityI[celli] = 1.0/max(condTable_m.value(Ti), SMALL);
    }

    kappa.correctBoundaryConditions();
    cp.correctBoundaryConditions();
    electrical_resistivity.correctBoundaryConditions();

    if (mthd.valid())
    {
        // MHD sees the alpha blended mixture like before, now at the cell T
        scalarField& condI = electrical_conductivity;
        scalarField& permI = magnetic_permeability;
        forAll(condI, celli)
        {
            const scalar a = max(min(alpha1I[celli], 1.0), 0.0);
            const scalar Ti = TI[celli];

            condI[celli] =
                a*condTable_m.value(Ti)
              + (1.0 - a)*condTable_gPtr->value(Ti);

            permI[celli] =
                a*magpermTable_mPtr->value(Ti)
              + (1.0 - a)*magpermTable_gPtr->value(Ti);
        }

        electrical_conductivity.correctBoundaryConditions();
        magnetic_permeability.correctBoundaryConditions();
    }
}
