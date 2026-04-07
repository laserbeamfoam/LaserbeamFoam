forAll(mesh.C(), celli)
{
    State[celli] = -2.0;
    transitionState[celli] = 0.0;

    if (condensateFiltered[celli] > 0.05)
    {
        const scalar Tsolidus = TSolidus[celli];
        const scalar Tliquidus = TLiquidus[celli];
        const scalar TvapCell = Tvap[celli];
        const scalar Tcell = T[celli];
        const scalar dTdt = heatingRate[celli];

        if (Tcell < Tsolidus)
        {
            State[celli] = 2.0;
        }
        else if (Tcell < Tliquidus)
        {
            State[celli] = 1.0;

            if (dTdt > 0.0)
            {
                transitionState[celli] = -1.0;
            }
            else if (dTdt < 0.0)
            {
                transitionState[celli] = 1.0;
            }
        }
        else if (Tcell < TvapCell)
        {
            State[celli] = 0.0;
        }
        else
        {
            State[celli] = -1.0;
        }
    }
}

State.correctBoundaryConditions();
transitionState.correctBoundaryConditions();
