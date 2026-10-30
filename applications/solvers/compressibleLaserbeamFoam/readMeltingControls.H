IOdictionary solutionDict
(
    IOobject
    (
        "fvSolution",
        runTime.system(), //folder where the fvSolution file is located
        runTime,
        IOobject::MUST_READ
    )
);



dictionary& meltingDict = solutionDict.subDict("MELTING");

int minTCorr(readLabel(meltingDict.lookup("minTempCorrector")));
int maxTCorr(readLabel(meltingDict.lookup("maxTempCorrector")));

dimensionedScalar epsilonTol(readScalar(meltingDict.lookup("epsilonTolerance")));
dimensionedScalar epsilonRel(readScalar(meltingDict.lookup("epsilonRelaxation")));

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