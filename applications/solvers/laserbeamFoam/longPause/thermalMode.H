// Thermal mode: only the energy equation (TEqn.H) is solved, with the
// interface frozen and the fluid at rest. It is entered when all the metal is
// solid and the lasers are off (see checkThermalEntry.H): with no flow and no
// laser, the convection and laser terms of TEqn.H are zero, and the liquid
// fraction does not change, so its latent heat term is zero too. The time
// step is limited by the temperature change per step and by the time at
// which the next laser is turned on. The function objects are executed and
// the fields written as in the full mode.
{
    Info<< "\nlongPause: thermal mode from t = " << runTime.timeName()
        << ", with the interface frozen and the fluid at rest\n" << endl;

    const scalar fullDeltaT = runTime.deltaTValue();

    // No laser and no flow
    laser.resetDeposition();
    U = dimensionedVector(U.dimensions(), Zero);
    U.correctBoundaryConditions();
    phi = dimensionedScalar(phi.dimensions(), Zero);
    rhoPhi = dimensionedScalar(rhoPhi.dimensions(), Zero);

    // The interface does not move in this mode: the properties that depend on
    // it are set once
    gradAlpha = fvc::grad(alpha1);
    gradAlpha.correctBoundaryConditions();

    #include "updateProps.H"

    scalar maxTemperatureChange = -1;
    bool resumeFull = false;

    while (!resumeFull && runTime.run())
    {
        #include "readControls.H"

        const scalar t = runTime.value();

        // Time step: initialDeltaT first, then grow by at most
        // maxDeltaTGrowth and limit the largest temperature change of the
        // previous step to maxTemperatureChange
        scalar deltaT = thermalInitialDeltaT;
        if (maxTemperatureChange >= 0)
        {
            deltaT = runTime.deltaTValue()*thermalMaxDeltaTGrowth;
            if (maxTemperatureChange > SMALL)
            {
                deltaT = min
                (
                    deltaT,
                    runTime.deltaTValue()
                   *thermalMaxTemperatureChange/maxTemperatureChange
                );
            }
        }
        deltaT = min(deltaT, thermalMaxDeltaT);

        // Equal steps to the end time, so that the last one ends on it
        const scalar remaining = runTime.endTime().value() - t;
        deltaT = remaining/std::ceil(remaining/deltaT);

        // Stop leadTime before the next laser is turned on
        const scalar tResume = laser.nextOnTime(t) - thermalLeadTime;
        if (t + deltaT >= tResume)
        {
            deltaT = max(tResume - t, SMALL);
            resumeFull = true;
        }

        // Without adjustment to the write times, which could lengthen the
        // step beyond the time at which the full mode must resume
        runTime.setDeltaT(deltaT, false);
        ++runTime;

        Info<< "Time = " << runTime.timeName() << " (thermal mode)" << nl
            << "deltaT = " << runTime.deltaTValue() << endl;

        #include "updateKappaCp.H"

        const scalarField T0(T.primitiveField());

        #include "TEqn.H"

        maxTemperatureChange = gMax(mag(T.primitiveField() - T0));

        Info<< "Largest temperature change = " << maxTemperatureChange
            << " K, min/max(T) = " << gMin(T) << ", " << gMax(T) << endl;

        runTime.write();

        // Write the final time also when it is not a write time (not at the
        // return to the full mode, which may be close to the end time)
        if
        (
            !resumeFull
         && runTime.value()
          > runTime.endTime().value() - 0.5*runTime.deltaTValue()
         && !runTime.writeTime()
        )
        {
            runTime.writeNow();
        }

        runTime.printExecutionTime(Info);
    }

    if (resumeFull)
    {
        // Continue with the time step of the last step in the full mode
        runTime.setDeltaT(fullDeltaT, false);

        Info<< "\nlongPause: full mode from t = " << runTime.timeName()
            << ", " << thermalLeadTime
            << " s before the next laser is turned on\n" << endl;
    }
    else
    {
        // The end time was reached in the thermal mode, and runTime.run()
        // has already executed and ended the function objects
        endReached = true;
    }
}
