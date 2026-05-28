{
    // Temporarily making U relative to mesh motion
    if (mesh.moving())
    {
        U -= fvc::reconstruct(mesh.phi());
    }

    // Updating alpha1
    #include "alphaSuSp.H"
    advector.advect(Sp, Su);

    // Making U absolute again after advection step
    if (mesh.moving())
    {
        U += fvc::reconstruct(mesh.phi());
    }

    #include "rhofs.H"
    rhoPhi = advector.getRhoPhi(rho1f, rho2f);

    alpha2 = 1.0 - alpha1;
    mixture.correct();
}

scalar domainFraction = 0;
if (porosityEnabled)
{
    const volScalarField& porosity = tporosity.cref();
    rhoPhi *= scalar(1)/fvc::interpolate(porosity);
    domainFraction = alpha1.weightedAverage(mesh.Vsc()*porosity).value();
}
else
{
    domainFraction = alpha1.weightedAverage(mesh.Vsc()).value();
}

Info<< "Phase-1 volume fraction = "
    << domainFraction
    << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
    << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
    << endl;
