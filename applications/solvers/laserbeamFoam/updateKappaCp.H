{
    const scalarField& alpha1I = alpha1;
    const scalarField& TI = T;
    scalarField& kappaI = kappa;
    scalarField& cpI = cp;
    forAll(kappaI, celli)
    {
        kappaI[celli] =
            (
                max(min(alpha1I[celli], 1.0), 0.0)
               *(polykappa_m.value(TI[celli]))
            )
          + (
                (1.0 - max(min(alpha1I[celli], 1.0), 0.0))
               *polykappa_g.value(TI[celli])
            );

        // The specific heat is per unit mass: mixed with the mass fraction of
        // the metal, so that rho*cp is the heat capacity of the two phases in
        // the cell
        const scalar a = max(min(alpha1I[celli], 1.0), 0.0);
        const scalar massFraction =
            a*rho1.value()/(a*rho1.value() + (1.0 - a)*rho2.value());

        cpI[celli] =
            massFraction*polycp_m.value(TI[celli])
          + (1.0 - massFraction)*polycp_g.value(TI[celli]);
    }

    kappa.correctBoundaryConditions();
    cp.correctBoundaryConditions();
}
