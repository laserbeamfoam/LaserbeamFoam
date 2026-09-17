{
    scalarField& eResI = electrical_resistivity;
    const scalarField& TI = T;
    const scalarField& alpha1I = alpha1;
    
    forAll(eResI, celli)
    {
    
    	scalar a = max(min(alpha1I[celli], 1.0), 0.0);

        eResI[celli] = a*electricalResistivity_m->value(TI[celli]);
            
    }    

    electrical_resistivity.correctBoundaryConditions();
    
    Info<< "electrical_resistivity min/max: "
    	<< min(electrical_resistivity).value() << " "
    	<< max(electrical_resistivity).value() << nl;
    
    
}


// end


