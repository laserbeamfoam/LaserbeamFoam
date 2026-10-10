{
    TRHS *= 0.0;

    rhoCp = rho*cp;
    rhoCp.correctBoundaryConditions();

    rhophicp = fvc::interpolate(cp)*rhoPhi;

    em = (epsilon1*emL) + ((1.0 - epsilon1)*emS);

    TRHS =
        LatentHeat*((fvc::ddt(rho, epsilon1) + fvc::div(rhoPhi, epsilon1)));
    TRHS.correctBoundaryConditions();

    Qv =
        0.82*LatentHeatVap*Mm*p0
       *Foam::exp(LatentHeatVap*Mm*((T - Tvap)/(R*T*Tvap)))
       /Foam::pow(2*M_PI*Mm*R*T, 0.5);
    Qv.correctBoundaryConditions();

    const volScalarField Qrad_expl(SB*em*(Foam::pow(TRef, 4) + 3*Foam::pow(T, 4)));
    const volScalarField Qrad_inpl(-(4*SB*em*Foam::pow(T, 3)));

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

    const volScalarField gradAlphaDamper(mag(gradAlpha)*thermalDamper);

    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T)
      + fvm::div(rhophicp, T)
      - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhophicp), T)
      - fvm::laplacian(kappa, T)
      - fvm::Sp(Qrad_inpl*gradAlphaDamper, T)
     ==
        laser.deposition()
    //   - Qv*mag(gradAlpha)*thermalDamper
      + (Qrad_expl - Qv)*gradAlphaDamper
      - TRHS
    );

    TEqn.relax();
    TEqn.solve();

    T.correctBoundaryConditions();

    gradT = fvc::grad(T);
    gradT.correctBoundaryConditions();
}
