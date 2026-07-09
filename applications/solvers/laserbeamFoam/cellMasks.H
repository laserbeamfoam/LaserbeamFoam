forAll(mesh.C(), celli)
{
    State[celli] = -2.0;  // gas
    transitionState[celli] = 0.0; // 0.0 = no transition, -1.0 = heating, 1.0 = cooling

    if (alpha1[celli] > 0.05)
    {
        const scalar Tsolidus = TSolidus[celli];
        const scalar Tliquidus = TLiquidus[celli];
        const scalar Tcell = T[celli];
        const scalar dTdt = heatingRate[celli];

        if (Tcell < Tsolidus)
        {
            State[celli] = 2.0;
        }
        else if (Tcell < Tliquidus)
        {
            State[celli] = 1.0; // mush

            if (dTdt > 0.0)
            {
                transitionState[celli] = -1.0; // Mark as heating if temperature is increasing
            }
            else if (dTdt < 0.0)
            {
                transitionState[celli] = 1.0; // Mark as cooling if temperature is decreasing
            }
        }
        else if (Tcell < Tvap.value())
        {
            State[celli] = 0.0; // liquid
        }
        else
        {
            State[celli] = -1.0; // gas/vapour
        }
    }
}

State.correctBoundaryConditions();
transitionState.correctBoundaryConditions();
