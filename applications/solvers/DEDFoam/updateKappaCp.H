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

        cpI[celli] =
            (
                max(min(alpha1I[celli], 1.0), 0.0)
               *polycp_m.value(TI[celli])
            )
          + (
                (1.0 - max(min(alpha1I[celli], 1.0), 0.0))
               *polycp_g.value(TI[celli])
            );
    }

    kappa.correctBoundaryConditions();
    cp.correctBoundaryConditions();
}
