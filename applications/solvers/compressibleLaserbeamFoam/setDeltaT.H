if (adjustTimeStep)
{
    // --- user limits (optional; sensible defaults) -----------------------
    const scalar maxAcousticCo =
        runTime.controlDict().getOrDefault<scalar>("maxAcousticCo", 50.0);
    // vapour volume made per step over the cell's: with the implicit update a Courant number, not a stability limit
    const scalar maxPhaseChangeFraction =
        runTime.controlDict().getOrDefault<scalar>("maxPhaseChangeFraction", 0.25);
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


    // phase-change volume per step at the rate before the per-step volume cap (see solveAlphas)
    const scalar phaseChangeNum =
        gMax(mag(vDotDt.primitiveField())())*deltaT0;

    // Fickian interdiffusion of the species is explicit (solveAlphas): a cell stays bounded while
    // dt*sum_f(D*|Sf|/|d|)/V < 1 over its faces to other cells. Inactive when no pair diffuses
    const scalar maxDiffusionNo =
        runTime.controlDict().getOrDefault<scalar>("maxDiffusionNo", 0.5);

    scalar diffusionNum = 0.0;
    const scalar Dmax = mixture.maxInterfaceDiffusivity();
    if (Dmax > 0)
    {
        const surfaceScalarField& dc = mesh.surfaceInterpolation::deltaCoeffs();
        const surfaceScalarField& magSf = mesh.magSf();
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();

        scalarField sumAd(mesh.nCells(), 0.0);
        forAll(own, facei)
        {
            const scalar Ad = magSf[facei]*dc[facei];
            sumAd[own[facei]] += Ad;
            sumAd[nei[facei]] += Ad;
        }
        forAll(mesh.boundary(), patchi)
        {
            const fvPatch& fvp = mesh.boundary()[patchi];
            if (!fvp.coupled()) continue;

            const labelUList& fc = fvp.faceCells();
            const scalarField& magSfp = magSf.boundaryField()[patchi];
            const scalarField& dcp = dc.boundaryField()[patchi];
            forAll(fc, i)
            {
                sumAd[fc[i]] += magSfp[i]*dcp[i];
            }
        }

        diffusionNum = Dmax*deltaT0*gMax(sumAd/mesh.V().field());
    }

    // physical rates only: flow, interface, sound, phase change, diffusion. Each one shrinks with dt
    const scalar rCo       = maxCo                 /(CoNum          + smallT);
    const scalar rAlpha    = maxAlphaCo            /(alphaCoNum     + smallT);
    const scalar rAcoustic = maxAcousticCo         /(acousticCoNum  + smallT);
    const scalar rPhase    = maxPhaseChangeFraction/(phaseChangeNum + smallT);

    scalar limiter = min(min(rCo, rAlpha), min(rAcoustic, rPhase));

    if (Dmax > 0)
    {
        limiter = min(limiter, maxDiffusionNo/(diffusionNum + smallT));
    }

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

    // the smallest ratio is the one setting dt (< 1 means it is cutting it)
    Info<< "deltaT = " << runTime.deltaTValue()
        << "  (Co " << rCo << ", interface Co " << rAlpha
        << ", acoustic " << rAcoustic << ", phase change " << rPhase;
    if (Dmax > 0)
    {
        Info<< ", diffusion " << maxDiffusionNo/(diffusionNum + smallT);
    }
    Info<< ")" << endl;
}
