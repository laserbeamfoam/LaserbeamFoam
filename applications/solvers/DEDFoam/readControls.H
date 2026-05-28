IOdictionary solutionDict
(
    IOobject
    (
        "fvSolution",
        runTime.system(),
        runTime,
        IOobject::MUST_READ
    )
);

dictionary& meltingDict = solutionDict.subDict("MELTING");

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
