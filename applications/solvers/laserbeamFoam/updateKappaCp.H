{
    const scalarField& alpha1I = alpha1;
    const scalarField& massFractionI = massFraction;
    const scalarField& TI = T;
    scalarField& kappaI = kappa;
    scalarField& cpI = cp;
    forAll(kappaI, celli)
    {
        const scalar a = max(min(alpha1I[celli], 1.0), 0.0);

        kappaI[celli] =
            a*polykappa_m.value(TI[celli])
          + (1.0 - a)*polykappa_g.value(TI[celli]);

        // cp is specific (per unit mass), so it is mixed with the mass
        // fraction so that rho*cp is volume-additive
        cpI[celli] =
            massFractionI[celli]*polycp_m.value(TI[celli])
          + (1.0 - massFractionI[celli])*polycp_g.value(TI[celli]);
    }

    kappa.correctBoundaryConditions();
    cp.correctBoundaryConditions();
}
