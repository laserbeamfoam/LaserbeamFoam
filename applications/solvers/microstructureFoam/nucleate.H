Info << "Nucleating cells: nucleate.H" << endl;

Tu = TLiquidus - T; // Calculate undercooling for convenience
Info << "Tu max: " << gMax(Tu) << "\tmin: " << gMin(Tu) << "\tavg: " << gAverage(Tu) << endl;
Info << "TLiquidus max: " << gMax(TLiquidus) << "\tmin: " << gMin(TLiquidus) << "\tavg: " << gAverage(TLiquidus) << endl;


//- DEACTIVATION LOOP
label sitesDeactivated(0);
//- Need to deactivate any sites are within solid
forAll(nucSite,i){
    // Deactivate any sites that are now within solid
    if (mag(nucSite[i]-1.0) < SMALL) // currently active site. Do I even need to do this?
    	// Why not just apply to all sites?
    {
        // GAS
        if ( alpha1[i] < 0.5){ // alpha1 = 1 for metal, 0 for gas
            nucSite[i] = 0.0; // Turn off this site
            sitesDeactivated++;
        }

        // SOLID
        else if (epsilon1mask[i] < 0.5) // NOT SOLID. 
        // epsilon1mask solid = 0, everything else = 1)
        {
            nucSite[i] = 0.0; // Turn off this site
            sitesDeactivated++;
        }
    }
}

//- ACTIVATION LOOP
label sitesActivated(0);
// Only sites within liquid metal should be activated.
forAll(nucSite,i){
    if (mag(nucSite[i]) < SMALL) // Currently inactive site (i.e. == 0)
    {
        if (
            (alpha1[i] > 0.5) // indicates metal
        &&  (epsilon1mask[i] > 0.5) // indicates not solid (must be liquid))
        )
        {
            nucSite[i] = 1.0;
            sitesActivated++;
        }
    }
}

// Need to check if it is nucleating somewhere an existing interface is
volScalarField globalSumNi = PopBal[0];

for (int i = 1; i < PopBal.size(); ++i)
{
    globalSumNi += PopBal[i];
}

// Use a dynamic list on each processor to record cells that should be nucleated.
DynamicList<label> cellsToBeNuced (0);
label numSitesToBeNuced(0);

forAll(nucSite,i){
    if (nucSite[i] > 0.5){ // Site is currently active
        if (  	(nucTu[i] >= 0) // Double check that it has a nucleation undercooling assigned, also, this condition is fastest to check.  
        	&&	(Tu[i] >= nucTu[i]) // Check if undercooling has been reached.
            &&  (globalSumNi[i] < SMALL) // Check if this location is free on all fields.
            )
        {
            // Cell is to be nucleated
            cellsToBeNuced.append(i);
            numSitesToBeNuced++;
            nucSite[i] = 0; // 0 = deactived nucleation site
        }
    }
}
Pout << "Proc " << Pstream::myProcNo() << "\tnNuc: " << numSitesToBeNuced << "\tnDeact: " << sitesDeactivated <<  "\tnAct: " << sitesActivated <<endl;

// int dummy;
// Pstream::scatter(dummy);

// abort();
// Now we have a list of cells that should be nucleated on each processor.
// Need to cycle through processors and nucleate all cells on the processor,
// updating maxNiIndex each time. When moving to next processor, pass on maxNiIndex

// The problem is write here:
//  In the reduce operation probably.

// The problem might be occuring when a cell on more than one proc nucleates.
// The problem is then confused at to which field it should be nucleating as
// the index of the latest assigned field is not sync'd

// Solution:
// 		Cycle through processors, one at a time when nucleating.

// Need a for loop to cycle through processors
// When nucleating 

reduce(maxNiIndex, maxOp<label>());

label nProcs = Pstream:: nProcs();
label myProc = Pstream::myProcNo();

// We want 
for(int i = 0; i < nProcs; i++){
    if (i == myProc)
    {
        printf("Proc %d wants to nucleate %d cells\n", myProc, cellsToBeNuced.size());
        forAll(cellsToBeNuced, j){
            label nucCellI = cellsToBeNuced[j];
            maxNiIndex++; // Increment, so that it nucleates it in the next available field
            PopBal[maxNiIndex][nucCellI] = 1.0;
            xi[nucCellI] = 1.0;
            printf("\nCell %d was nucleated by proc %d on field n.%d\n\n", nucCellI, myProc, maxNiIndex);
        }
    }
    // Need to sync this new maxNiIndex with the other processes
    reduce(maxNiIndex, maxOp<label>());   
}

// forAll(cellsToBeNuced, j){
//     label nucCellI = cellsToBeNuced[j];
//     maxNiIndex++; // Increment, so that it nucleates it in the next available field
//     PopBal[maxNiIndex][nucCellI] = 1.0;
//     xi[nucCellI] = 1.0;
//     // reduce(maxNiIndex, maxOp<label>());
// }


Info << "end of nucleation nProcs: " << Pstream::nProcs() << endl;