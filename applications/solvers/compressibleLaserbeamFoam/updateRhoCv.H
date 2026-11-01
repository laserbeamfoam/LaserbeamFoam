// Mass-weighted mixture heat capacity for the energy equation (TEqn.H).
//
//  rhoCv     = sum_i      alpha_i*rho_i*Cv_i   [J/m^3/K]  volume-additive
//  rhoCvCond = sum_i,cond alpha_i*rho_i*Cv_i   [J/m^3/K]  condensed only
//  CvMass    = rhoCv/rho = sum_i Y_i*Cv_i      [J/kg/K]   Y_i = alpha_i rho_i/rho
//
// Cv is a per-unit-mass property, so the only mixing rule that makes
// rho*Cv volume-additive is the MASS-fraction one. The old mixture.rCv()
// (= sum_i alpha_i/Cv_i, volume weighted) is what compressibleInterFoam
// uses; it is only exact in pure-phase cells. In an interface cell it
// overcounts every source (laser, Lvap*r, Joule, work) by the factor
// CvMass*sum_i(alpha_i/Cv_i), ~1.5-2 for metal/metal-vapour cells.
//
// Called at the top of every outer corrector (TEqn.H) because rho_i(p,T)
// and Cv_i(T) move inside the step, and once in createFields.H to seed the
// oldTime storage. Uses the solver field rho (= mixture.rho() at the last
// update) so CvMass*rho == rhoCv exactly, cell by cell.
{
    rhoCv == dimensionedScalar(rhoCv.dimensions(), Zero);
    rhoCvCond == dimensionedScalar(rhoCvCond.dimensions(), Zero);

    for (const phaseModel& ph : mixture.phases())
    {
        const volScalarField rhoCvPh(ph*ph.thermo().rho()*ph.thermo().Cv());

        rhoCv += rhoCvPh;

        if (ph.isCondensed())
        {
            rhoCvCond += rhoCvPh;
        }
    }

    rhoCv.correctBoundaryConditions();
    rhoCvCond.correctBoundaryConditions();

    CvMass = rhoCv/max(rho, dimensionedScalar("rhoSmallCv", dimDensity, SMALL));
    CvMass.correctBoundaryConditions();
}
