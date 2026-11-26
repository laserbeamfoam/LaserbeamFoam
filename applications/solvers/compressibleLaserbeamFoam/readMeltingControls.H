// MELTING controls, from the mesh's own fvSolution: no file is read here, and an edited fvSolution is
// picked up by the run-time re-read as for any other solver control
const dictionary& meltingDict = mesh.solutionDict().subDict("MELTING");

const label minTCorr(meltingDict.get<label>("minTempCorrector"));
const label maxTCorr(meltingDict.get<label>("maxTempCorrector"));

const dimensionedScalar epsilonTol(meltingDict.get<scalar>("epsilonTolerance"));
const dimensionedScalar epsilonRel(meltingDict.get<scalar>("epsilonRelaxation"));

// Treatment of the latent heat of fusion in the melt-fraction correction loop
// (same keyword and meaning as laserbeamFoam PR 159):
//   apparentCapacity (default)

const word latentHeatMethod
(
    meltingDict.getOrDefault<word>("latentHeatMethod", "apparentCapacity")
);

if
(
    latentHeatMethod != "sourceBased"
 && latentHeatMethod != "apparentCapacity"
)
{
    FatalIOErrorInFunction(meltingDict)
        << "latentHeatMethod must be sourceBased or apparentCapacity, not "
        << latentHeatMethod << exit(FatalIOError);
}
