// Decide whether the thermal mode can be entered at the current time:
// all lasers off with no laser turned on for at least minPause, no liquid
// metal, metal below the solidus by solidusMargin, and metal at rest.
// Without a heat source, conduction cannot raise the temperature above its
// current maximum, so the metal cannot melt during the thermal mode.

bool enterThermal = false;
{
    const scalar t = runTime.value();

    if (laser.allOff(t) && laser.nextOnTime(t) - t > thermalMinPause)
    {
        // Cells whose solidus and liquidus are those of the metal
        // (see updateProps.H)
        scalar maxLiquidFraction = 0;
        scalar maxMetalT = 0;
        scalar maxMetalU = 0;

        const scalarField& alpha1I = alpha1;
        const scalarField& epsilon1I = epsilon1;
        const scalarField& TI = T;
        const vectorField& UI = U;

        forAll(alpha1I, celli)
        {
            if (alpha1I[celli] > 0.05)
            {
                maxLiquidFraction = max(maxLiquidFraction, epsilon1I[celli]);
                maxMetalT = max(maxMetalT, TI[celli]);
            }

            if (alpha1I[celli] > 0.5)
            {
                maxMetalU = max(maxMetalU, mag(UI[celli]));
            }
        }

        reduce(maxLiquidFraction, maxOp<scalar>());
        reduce(maxMetalT, maxOp<scalar>());
        reduce(maxMetalU, maxOp<scalar>());

        enterThermal =
            maxLiquidFraction <= thermalMaxLiquidFraction
         && maxMetalT < Tsolidus1.value() - thermalSolidusMargin
         && maxMetalU <= thermalMaxMetalVelocity;

        Info<< "longPause check: max metal liquid fraction = "
            << maxLiquidFraction << ", max metal T = " << maxMetalT
            << ", max metal |U| = " << maxMetalU << endl;
    }
}
