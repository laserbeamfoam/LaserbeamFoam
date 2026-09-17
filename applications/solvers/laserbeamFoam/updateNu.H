{
    scalarField& nuI = nuTableField;
    
    const scalarField& TI = T;
    const scalarField& alpha1I = alpha1;
    
    forAll(nuI, celli)
    {
     	scalar a = max(min(alpha1I[celli], 1.0), 0.0);
		scalar Tc = TI[celli];
		
        nuI[celli] = a*nu_m->value(Tc) + (1.0 - a)*nu_g->value(Tc);
            
    }    

    nuTableField.correctBoundaryConditions();
    
    Info<< "nuTableField min/max: "
    	<< min(electrical_resistivity).value() << " "
    	<< max(electrical_resistivity).value() << nl;

}


// end

