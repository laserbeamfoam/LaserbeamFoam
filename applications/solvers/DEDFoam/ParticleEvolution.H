// ============================================================================
// ParticleEvolution.H
//
// Evolves the kinematic particle cloud and accounts for particles melting
// above particle_Tmelt. Each parcel's mass is tracked exclusively through
// its diameter p.d(), which the lagrangian framework already serialises for
// restart. No external bookkeeping is needed.
//
// Melt model (exponential decay of current mass):
//     dm/dt = -m / time_melt
// integrated forward Euler. time_melt is the 1/e timescale (~37% mass left
// after time_melt; ~5% after 3*time_melt; effectively gone after ~5*time_melt).
//
// Stick condition: parcels in cells with alpha1 > alphaStickThreshold are
// re-velocitied to the local fluid velocity, so they entrain with the bath
// instead of ballistically passing through the molten pool into the solid
// substrate below.
//
// Sources delivered each step:
//   - alphaMeltSource           [1/s]      volume-fraction addition rate
//                                          (applied via MULES Su and a
//                                           div(U) source in pEqn)
//   - particleEnthalpySource    [W/m^3]    energy sink (RHS of TEqn)
//   - massSource                [kg/m^3/s] diagnostic mass-addition rate
//
// Restart: deleting the Map<scalar> made this naturally restart-safe.
// Diameters of any partially-melted particles are read from the lagrangian
// time directory; melting continues seamlessly from there.
// ============================================================================
{
    const scalar dt = runTime.deltaTValue();

    // Tunables — could be promoted to kinematicCloudProperties later.
    const scalar alphaStickThreshold = 0.5;

    // Kill parcels whose diameter falls below this threshold. Below
    // roughly the smallest cell size the parcel is no longer
    // meaningful, and continuing to track it wastes time. The
    // remaining mass is dumped into the bath in the killing step so
    // mass is conserved across the kill event.
    const scalar particleDeathDiameter = 5.0e-6;   // 5 μm
    const scalar particleDeathMass =
        constant::mathematical::pi/6.0
       *Foam::pow(particleDeathDiameter, 3)*rho1.value();

    // Reset source fields
    alphaMeltSource        *= 0.0;
    particleEnthalpySource *= 0.0;
    massSource             *= 0.0;

    label nMelted          = 0;
    label nPartiallyMelted = 0;
    label nStuck           = 0;

    Info<< "Evolving " << parcels.name() << endl;
    parcels.evolve();

    // ------------------------------------------------------------------
    // DIAGNOSTIC: where are the particles, and is any cell hot enough?
    // Remove or comment out once the melt source is producing values
    // you expect.
    // ------------------------------------------------------------------
    {
        label nActive       = 0;
        label nOutOfMesh    = 0;
        label nInHotCell    = 0;     // T[celli] > Tmelt
        label nInMetalCell  = 0;     // alpha1[celli] > 0.5
        label nHotAndInMetal= 0;
        scalar yMin =  GREAT, yMax = -GREAT;
        scalar Tmin =  GREAT, Tmax = -GREAT;  // T at occupied cells
        scalar alphaMin =  GREAT, alphaMax = -GREAT;  // alpha at occupied cells
        scalar dMin =  GREAT, dMax = -GREAT;  // diameter range of active parcels

        // State-persistence observables. Spawn comes from a Rosin-Rammler
        // truncated to [60, 100] um (set in kinematicCloudProperties).
        // Anything below dShrunkThreshold has undergone irreversible
        // melting at some point in its history. Counting these directly
        // proves diameters persist across timesteps.
        const scalar dShrunkThreshold = 59.0e-6;   // just under spawn min
        label nShrunk     = 0;
        scalar dSum       = 0.0;        // for mean diameter

        for (const basicKinematicParcel& p : parcels)
        {
            if (!p.active()) continue;
            nActive++;

            const label celli = p.cell();
            if (celli < 0 || celli >= mesh.nCells())
            {
                nOutOfMesh++;
                continue;
            }

            const scalar y = p.position().y();
            yMin = min(yMin, y);
            yMax = max(yMax, y);

            const scalar Tc = T[celli];
            const scalar ac = alpha1[celli];
            Tmin = min(Tmin, Tc);
            Tmax = max(Tmax, Tc);
            alphaMin = min(alphaMin, ac);
            alphaMax = max(alphaMax, ac);

            // Also track diameter range so we can see whether melting
            // is reflected in the d field across timesteps.
            const scalar dp = p.d();
            dMin = min(dMin, dp);
            dMax = max(dMax, dp);
            dSum += dp;
            if (dp < dShrunkThreshold) nShrunk++;

            const bool hot   = (Tc > particle_Tmelt.value());
            const bool metal = (ac > 0.5);
            if (hot)   nInHotCell++;
            if (metal) nInMetalCell++;
            if (hot && metal) nHotAndInMetal++;
        }

        reduce(nActive,         sumOp<label>());
        reduce(nOutOfMesh,      sumOp<label>());
        reduce(nInHotCell,      sumOp<label>());
        reduce(nInMetalCell,    sumOp<label>());
        reduce(nHotAndInMetal,  sumOp<label>());
        reduce(yMin,            minOp<scalar>());
        reduce(yMax,            maxOp<scalar>());
        reduce(Tmin,            minOp<scalar>());
        reduce(Tmax,            maxOp<scalar>());
        reduce(alphaMin,        minOp<scalar>());
        reduce(alphaMax,        maxOp<scalar>());
        reduce(dMin,            minOp<scalar>());
        reduce(dMax,            maxOp<scalar>());
        reduce(nShrunk,         sumOp<label>());
        reduce(dSum,            sumOp<scalar>());

        const scalar dMean =
            (nActive > 0) ? dSum/scalar(nActive) : 0.0;

        Info<< "ParticleDiag: active="     << nActive
            << ", out-of-mesh="            << nOutOfMesh
            << ", in-hot-cell="            << nInHotCell
            << ", in-metal-cell="          << nInMetalCell
            << ", in-both="                << nHotAndInMetal
            << "; y range=["               << yMin << "," << yMax << "]"
            << "; T at parcels=["          << Tmin << "," << Tmax << "]"
            << "; alpha at parcels=["      << alphaMin << "," << alphaMax << "]"
            << "; d range=["               << dMin*1e6 << "," << dMax*1e6 << "] um"
            << "; d mean="                 << dMean*1e6 << " um"
            << "; n with d<"               << dShrunkThreshold*1e6 << "um="
                                           << nShrunk
            << "; Tmelt="                  << particle_Tmelt.value()
            << endl;
    }

    // ------------------------------------------------------------------
    // Stick condition: any parcel that ended up inside the metal phase
    // gets its velocity replaced by the local mixture velocity so it
    // entrains with the bath rather than coasting through. This is
    // applied AFTER evolve(), so it affects the NEXT step's motion.
    // ------------------------------------------------------------------
    for (basicKinematicParcel& p : parcels)
    {
        if (!p.active()) continue;

        const label celli = p.cell();
        if (celli < 0 || celli >= mesh.nCells()) continue;

        if (alpha1[celli] > alphaStickThreshold)
        {
            p.U() = U[celli];
            nStuck++;
        }
    }

    // ------------------------------------------------------------------
    // First pass: aggregate desired melt per cell.
    //   dm = m * dt / time_melt   (current mass, not spawn mass)
    // ------------------------------------------------------------------
    scalarField proposedAlphaAddition (mesh.nCells(), 0.0);
    scalarField proposedMassAddition  (mesh.nCells(), 0.0);
    scalarField proposedEnthalpySink  (mesh.nCells(), 0.0);

    for (basicKinematicParcel& p : parcels)
    {
        if (!p.active()) continue;

        const label celli = p.cell();
        if (celli < 0 || celli >= mesh.nCells()) continue;

        const scalar tempP = T[celli];
        if (tempP <= particle_Tmelt.value()) continue;

        const scalar currentMass = p.mass();
        if (currentMass < SMALL || dt < SMALL) continue;

        // Forward-Euler step on dm/dt = -m/tau, then check whether the
        // result would put us below the death threshold. If so, dump
        // all remaining mass into the cell this step. Otherwise just
        // take the exponential step.
        const scalar exponentialMelt =
            min(currentMass, currentMass*dt/time_melt.value());
        const scalar massAfter = currentMass - exponentialMelt;

        const scalar meltedMass =
            (massAfter < particleDeathMass) ? currentMass : exponentialMelt;

        const scalar volumeAsMetal = meltedMass/rho1.value();
        const scalar cellVolume    = mesh.V()[celli];
        const scalar deltaAlpha    = volumeAsMetal/cellVolume;

        proposedAlphaAddition[celli] += deltaAlpha;
        proposedMassAddition [celli] += meltedMass;

        const scalar sensibleHeat =
            meltedMass*CpParticle.value()*(tempP - Tinject.value());
        const scalar latentHeat   = meltedMass*LfParticle.value();

        proposedEnthalpySink[celli] += (sensibleHeat + latentHeat);
    }

    // ------------------------------------------------------------------
    // Second pass: write the source fields directly. No per-cell space
    // limiter — that was the v1 approach, before pEqn carried a
    // div(U) = alphaMeltSource term. With the dilation source in place,
    // a fully-saturated cell (alpha1 = 1) accommodates the added metal
    // by velocity outflow, not by exceeding alpha1 = 1 locally. MULES
    // enforces 0 <= alpha1 <= 1 anyway, so the upper bound is hard.
    //
    // Excessively large proposals would be silently clipped by MULES
    // (small mass leak). The diagnostic at the bottom of this file
    // prints a warning if max predicted alpha after the step exceeds 1.
    // ------------------------------------------------------------------
    forAll(proposedAlphaAddition, celli)
    {
        if (proposedAlphaAddition[celli] <= SMALL) continue;

        const scalar cellVolume = mesh.V()[celli];

        alphaMeltSource[celli] =
            proposedAlphaAddition[celli]/dt;

        massSource[celli] =
            proposedMassAddition[celli]/(cellVolume*dt);

        particleEnthalpySource[celli] =
          - proposedEnthalpySink[celli]/(cellVolume*dt);
    }

    // ------------------------------------------------------------------
    // Third pass: shrink each parcel's diameter by the full desired
    // melt amount and kill parcels whose remaining mass is below SMALL.
    //
    // Uses a range-based for so that p is unambiguously a non-const
    // reference to the parcel stored in the cloud — some OpenFOAM
    // versions had subtleties with forAllIter that could silently
    // produce a temporary copy here, which would make the diameter
    // update vanish.
    // ------------------------------------------------------------------
    static label thirdPassDiagCount = 0;
    const label thirdPassDiagLimit = 5;

    for (basicKinematicParcel& p : parcels)
    {
        if (!p.active()) continue;

        const label celli = p.cell();
        if (celli < 0 || celli >= mesh.nCells()) continue;

        const scalar tempP = T[celli];
        if (tempP <= particle_Tmelt.value()) continue;

        const scalar currentMass = p.mass();
        if (currentMass < SMALL || dt < SMALL) continue;

        const scalar exponentialMelt =
            min(currentMass, currentMass*dt/time_melt.value());
        const scalar massAfter = currentMass - exponentialMelt;

        const scalar dBefore = p.d();

        if (massAfter < particleDeathMass)
        {
            // The exponential step would put the parcel below the
            // death threshold; kill it (mass already dumped to the
            // alpha source by the first pass).
            p.active(false);
            nMelted++;
        }
        else
        {
            const scalar rhoParticle = p.rho();
            const scalar newDiameter =
                Foam::pow
                (
                    6.0*massAfter/(rhoParticle*constant::mathematical::pi),
                    1.0/3.0
                );
            p.d() = newDiameter;
            nPartiallyMelted++;
        }

        if (thirdPassDiagCount < thirdPassDiagLimit)
        {
            Info<< "  meltDiag: parcel origId=" << p.origId()
                << " origProc=" << p.origProc()
                << " currentMass=" << currentMass
                << " exponentialMelt=" << exponentialMelt
                << " massAfter=" << massAfter
                << " dBefore=" << dBefore
                << " dAfter=" << p.d()
                << " active=" << label(p.active())
                << endl;
            thirdPassDiagCount++;
        }
    }

    // Parallel reductions for reporting
    reduce(nMelted,          sumOp<label>());
    reduce(nPartiallyMelted, sumOp<label>());
    reduce(nStuck,           sumOp<label>());

    if (nMelted > 0 || nPartiallyMelted > 0 || nStuck > 0)
    {
        const scalar totalMassAdded =
            gSum(massSource.primitiveField()*mesh.V())*dt;

        Info<< "Particles: stuck-to-bath "      << nStuck
            << ", fully melted "                << nMelted
            << ", partially melted "            << nPartiallyMelted
            << ", total mass added "            << totalMassAdded*1e3 << " g"
            << endl;
    }

    // Diagnostic sanity check
    {
        const scalar maxAlphaSource = gMax(alphaMeltSource);
        if (maxAlphaSource > SMALL)
        {
            const scalar maxPredictedAlpha =
                gMax
                (
                    alpha1.primitiveField()
                  + alphaMeltSource.primitiveField()*dt
                );

            Info<< "Max alphaMeltSource: "       << maxAlphaSource    << " 1/s, "
                << "max predicted alpha after: " << maxPredictedAlpha << endl;

            if (maxPredictedAlpha > 1.0)
            {
                WarningInFunction
                    << "Alpha may exceed 1.0! Max predicted = "
                    << maxPredictedAlpha << endl;
            }
        }
    }

    // ------------------------------------------------------------------
    // Cap excessive enthalpy sinks: a single time-step should not be
    // allowed to drop the local cell temperature by more than
    // maxAllowedTempDrop.
    // ------------------------------------------------------------------
    {
        const scalar maxAllowedTempDrop = 200.0;  // K
        forAll(particleEnthalpySource, celli)
        {
            if (particleEnthalpySource[celli] < -SMALL)
            {
                const scalar maxCoolingPower =
                    rhoCp[celli]*maxAllowedTempDrop/dt;

                if (mag(particleEnthalpySource[celli]) > maxCoolingPower)
                {
                    particleEnthalpySource[celli] = -maxCoolingPower;
                }
            }
        }
    }
    Info<< "min particleEnthalpySource "
        << gMin(particleEnthalpySource) << " W/m^3" << endl;
}