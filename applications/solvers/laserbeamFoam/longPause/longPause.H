// Enter the thermal mode if the conditions of checkThermalEntry.H hold, and
// return to the full mode leadTime before the next laser is turned on. It is
// checked at start-up, so that a run restarted within a pause continues in the
// thermal mode, and at the end of every time step of the full mode.
{
    #include "checkThermalEntry.H"

    if (enterThermal)
    {
        #include "thermalMode.H"
    }
}
