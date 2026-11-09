// Mass-weighted heat capacity for TEqn.H: rhoCv = sum_i c_i*Cv_i [J/m^3/K] (rhoCvCond: condensed only),
// CvMass = rhoCv/rhoM. Cv is per mass, so only mass weighting is volume-additive (compressibleInterFoam's
// volume-weighted rCv over-counts interface-cell sources 1.5-2x). Every outer corrector: Cv_i(T) moves
{
    const dimensionedScalar cZeroCv("cZeroCv", dimDensity, Zero);

    rhoCv == dimensionedScalar(rhoCv.dimensions(), Zero);
    rhoCvCond == dimensionedScalar(rhoCvCond.dimensions(), Zero);

    label phasei = 0;
    for (const phaseModel& ph : mixture.phases())
    {
        volScalarField cEff("cEff." + ph.name(), max(mixture.c()[phasei], cZeroCv));

        if (ledgerFloorEOS)
        {
            cEff = max(cEff, ph*ph.thermo().rho());   // same guard as rhoM
        }

        const volScalarField rhoCvPh(cEff*ph.thermo().Cv());

        rhoCv += rhoCvPh;

        if (ph.isCondensed())
        {
            rhoCvCond += rhoCvPh;
        }

        ++phasei;
    }

    // positive diagonal for fvm::ddt(rhoCv, T) in empty cells; 100 J/kg/K is below any Cv in the dicts
    rhoCv = max
    (
        rhoCv,
        dimensionedScalar("rhoCvMin", rhoCv.dimensions(), 100.0*rhoMMin.value())
    );

    rhoCv.correctBoundaryConditions();
    rhoCvCond.correctBoundaryConditions();

    CvMass = rhoCv/rhoM;
    CvMass.correctBoundaryConditions();
}
