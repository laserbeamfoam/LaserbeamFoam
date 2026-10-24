// Use the fvSolution dictionary held by the mesh, which is re-read when the
// file is modified, rather than reading the file at every time step. As for
// the other fvSolution controls, changes made to the file during the run are
// therefore only picked up when runTimeModifiable is on.
const dictionary& meltingDict = mesh.solution().subDict("MELTING");

const label minTCorr(readLabel(meltingDict.lookup("minTempCorrector")));

const label maxTCorr(readLabel(meltingDict.lookup("maxTempCorrector")));

const dimensionedScalar epsilonTol
(
    readScalar(meltingDict.lookup("epsilonTolerance"))
);

const dimensionedScalar epsilonRel
(
    readScalar(meltingDict.lookup("epsilonRelaxation"))
);

const bool damperSwitch =
    meltingDict.lookupOrDefault<bool>("damperSwitch", false);

// Treatment of the latent heat in the correction loop of the liquid fraction:
// sourceBased (the explicit source and the correction
// epsilon1 += epsilonRelaxation*cp/L*(T - Tcorr) of each correction) or
// apparentCapacity (the latent heat of the cells in the mushy range taken
// implicitly in T, a Newton linearisation; see TEqn.H)
const word latentHeatMethod
(
    meltingDict.getOrDefault<word>("latentHeatMethod", "sourceBased")
);
if (latentHeatMethod != "sourceBased" && latentHeatMethod != "apparentCapacity")
{
    FatalIOErrorInFunction(meltingDict)
        << "latentHeatMethod must be sourceBased or apparentCapacity, not "
        << latentHeatMethod << exit(FatalIOError);
}
