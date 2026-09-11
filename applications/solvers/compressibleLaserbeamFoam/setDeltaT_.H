if (adjustTimeStep)
{
    // --- user limits (optional; sensible defaults) -----------------------
    const scalar maxAcousticCo =
        runTime.controlDict().getOrDefault<scalar>("maxAcousticCo", 50.0);
    const scalar maxPhaseChangeFraction =
        runTime.controlDict().getOrDefault<scalar>("maxPhaseChangeFraction", 0.05);
    const scalar maxDeltaTK =
        runTime.controlDict().getOrDefault<scalar>("maxDeltaTK", 50.0); // [K]
    const scalar dampFactor =
        runTime.controlDict().getOrDefault<scalar>("deltaTDampFactor", 1.15);

    const scalar deltaT0 = runTime.deltaTValue();
    const scalar smallT  = SMALL;

    scalar acousticCoNum = 0.0;
    {
        const volScalarField pSafe(max(p, dimensionedScalar(dimPressure, 1e4)));
const volScalarField c
(
    sqrt(max(mixture.gamma()*pSafe/max(mixture.rho(), dimensionedScalar(dimDensity, SMALL)),
             dimensionedScalar(sqr(dimVelocity), Zero)))
);

        surfaceScalarField cf(fvc::interpolate(c));
        acousticCoNum =
            0.5*gMax
            (
                (mesh.surfaceInterpolation::deltaCoeffs()*cf)().primitiveField()
            )*deltaT0;
    }


    scalar phaseChangeNum = gMax(mag(vDot.primitiveField()))*deltaT0;

    scalar deltaTKnum = 0.0;
    {
        const volScalarField gradPbyRho
        (
            mag(fvc::grad(p))/max(mixture.rho(), dimensionedScalar(dimDensity, smallT))
        );
        // ΔT_K per unit dt:  |U| * (|∇p|/rho) * rCv   [K/s]
        const volScalarField dTKdt
        (
            mag(U)*gradPbyRho*mixture.rCv()
        );
        deltaTKnum = gMax(dTKdt.primitiveField())*deltaT0;
    }


    const scalar rCo        = maxCo            /(CoNum         + smallT);
    const scalar rAlpha     = maxAlphaCo       /(alphaCoNum    + smallT);
    const scalar rAcoustic  = maxAcousticCo    /(acousticCoNum + smallT);
    const scalar rPhase     = maxPhaseChangeFraction/(phaseChangeNum + smallT);
    const scalar rDeltaTK   = maxDeltaTK       /(deltaTKnum    + smallT);

    const scalar limiter = min(min(min(rCo, rAlpha), min(rAcoustic, rPhase)), rDeltaTK);

    scalar deltaTFact = min(min(limiter, dampFactor*1.0), 1.2);

    runTime.setDeltaT
    (
        min(deltaTFact*deltaT0, maxDeltaT)
    );

    const scalar minDeltaT = runTime.controlDict().getOrDefault<scalar>("minDeltaT", 1e-12);
if (runTime.deltaTValue() < minDeltaT)
{
    runTime.setDeltaT(minDeltaT);
    Info<< "deltaT floored to minDeltaT = " << minDeltaT << endl;
}



    Info<< "deltaT = " << runTime.deltaTValue() << endl;
}