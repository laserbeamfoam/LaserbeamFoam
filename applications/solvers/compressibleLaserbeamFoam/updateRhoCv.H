// Mass-weighted mixture heat capacity for the energy equation (TEqn.H).
//
//  rhoCv     = sum_i      c_i*Cv_i    [J/m^3/K]  volume-additive
//  rhoCvCond = sum_i,cond c_i*Cv_i    [J/m^3/K]  condensed only
//  CvMass    = rhoCv/rhoM = sum_i Y_i*Cv_i  [J/kg/K]   Y_i = c_i/rhoM
//
// c_i = mixture.c()[i] is the ledger mass (alpha_i*rho_i only where the closure is exact)
//
// Cv is a per-unit-mass property, so the only mixing rule that makes
// rho*Cv volume-additive is the MASS-fraction one. The old mixture.rCv()
// (= sum_i alpha_i/Cv_i, volume weighted) is what compressibleInterFoam
// uses; it is only exact in pure-phase cells. In an interface cell it
// overcounts every source (laser, Lvap*r, Joule, work) by the factor
// CvMass*sum_i(alpha_i/Cv_i), ~1.5-2 for metal/metal-vapour cells.
//
// Called at the top of every outer corrector (TEqn.H) because Cv_i(T) moves
// inside the step, and once in createFields.H to seed the oldTime storage.
// Uses the solver field rhoM (frozen over the step, like c_i) so
// CvMass*rhoM == rhoCv exactly, cell by cell.
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
