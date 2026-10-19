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
