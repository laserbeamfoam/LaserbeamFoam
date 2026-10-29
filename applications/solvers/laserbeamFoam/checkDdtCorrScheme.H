// The ddtCorr flux is calculated once per outer corrector, before the pressure
// corrector loop, which is exact when it depends only on the fields of the
// previous time step. CoEuler and SLTS calculate a local time step from the
// current flux phi, which the pressure correctors update, so they are not
// supported; bounded is checked through the scheme it wraps. Checked in every
// outer corrector, as the schemes may be re-read.
{
    ITstream& ddtUScheme = mesh.ddtScheme("ddt(" + U.name() + ')');

    word ddtUSchemeName(ddtUScheme);
    if (ddtUSchemeName == "bounded")
    {
        ddtUScheme >> ddtUSchemeName;
    }
    ddtUScheme.rewind();

    if (ddtUSchemeName == "CoEuler" || ddtUSchemeName == "SLTS")
    {
        FatalErrorInFunction
            << "The ddt scheme " << ddtUSchemeName << " of ddt(" << U.name()
            << ") is not supported: it calculates a local time step from the "
            << "current flux, and laserbeamFoam calculates the ddtCorr flux "
            << "once per outer corrector, before the pressure correctors "
            << "update the flux. Use another ddt scheme, such as Euler."
            << exit(FatalError);
    }
}
