{


dimensionedScalar gamma_pf("gamma_pf", dimensionSet(0, 0, 0, 0, 0), readScalar(PhaseFieldProperties.lookup("gamma_pf")));

dimensionedScalar sigma_P("sigma_P", dimensionSet(1, 0, -2, 0, 0),readScalar(PhaseFieldProperties.lookup("sigma_P")));
dimensionedScalar sigma_G("sigma_G", dimensionSet(1, 0, -2, 0, 0),readScalar(PhaseFieldProperties.lookup("sigma_G")));
dimensionedScalar lp("lp", dimensionSet(0, 1, 0, 0, 0),readScalar(PhaseFieldProperties.lookup("lp")));
dimensionedScalar lg("lg", dimensionSet(0, 1, 0, 0, 0),readScalar(PhaseFieldProperties.lookup("lg")));

dimensionedScalar deltaf_g("deltaf_g", dimensionSet(0, 0, 0, 0, 0),readScalar(PhaseFieldProperties.lookup("deltaf_g")));
dimensionedScalar deltaf_p("deltaf_p", dimensionSet(0, 0, 0, 0, 0),readScalar(PhaseFieldProperties.lookup("deltaf_p")));

dimensionedScalar D0("D0", dimensionSet(-1, 1, 1, 0, 0),readScalar(PhaseFieldProperties.lookup("D0")));
dimensionedScalar Qg("Qg", dimensionSet(1, 2, -2, 0, -1),readScalar(PhaseFieldProperties.lookup("Qg")));


dimensionedScalar m_p = (3.0/4.0)*((sigma_P)/(deltaf_p*lp));
// dimensionedScalar m_p("m_p", dimensionSet(0, 0, -1, 0, 0), (3.0/4.0)*((sigma_P.value())/(deltaf_p.value()*lp.value())));

dimensionedScalar m_g= (3.0/4.0)*((sigma_G)/(deltaf_g*lg));
dimensionedScalar Kg = sigma_G*lg;
dimensionedScalar Kp = (3.0/4.0)*sigma_P*lp;


// dimensionedScalar Lg("Lg", dimensionSet(-1, 1, 1, 0, 0), readScalar(PhaseFieldProperties.lookup("Lg")));

dimensionedScalar gasconstant("gasconstant",dimensionSet(1, 2, -2, -1, -1),scalar(8.314));

// volScalarField 

if(Aniso_strength.value()>SMALL){
Lg = (GBs+Aniso_Mob)*(D0*Foam::exp((-Qg)/(gasconstant*T)))*1.0;
}
else{
Lg = (D0*Foam::exp((-Qg)/(gasconstant*T)))*1.0;
}
//Lg = (GBs+Aniso_Mob)*(D0*Foam::exp((-Qg)/(gasconstant*T)))*1.0;

// dimensionedScalar Lp("Lp", dimensionSet(-1, 1, 1, 0, 0), readScalar(PhaseFieldProperties.lookup("Lp")));

// dimensionedScalar 

if(Aniso_strength.value()>SMALL){
Lp = Aniso_Mob*(1.0*D0*Foam::exp((-Qg)/(gasconstant*T)));
}
else{
    Lp = (1.0*D0*Foam::exp((-Qg)/(gasconstant*T)));
}


Info << "Before sum phases " << endl;
Sum_Phases*=0.0;
GBs*=0.0;

/**/
for (label count = 0; count < N_Ori; count++) //find the number of non zero phases at each point in the domain
{
    volScalarField & ni = PopBal[count];
    PopBal_av[count] = fvc::average(ni); //takes the average from all cells sharing a face - this way the information about phases in adjacent cells is shared with the local cell
    PopBal_av[count].correctBoundaryConditions();
    Sum_parts[count]*=0.0;
}


// forAll(mesh.C(), cell)
// {
//     label active_phases[N_Ori] =
//     {
//         0
//     };
//     label n_active_phases = 0;



//     for (label count = 0; count < N_Ori; count++)
//     {
//         if (PopBal_av[count][cell] > SMALL)
//         {
//             active_phases[n_active_phases++] = count;
//         }

//     }


//     for (label count_i = 0; count_i < n_active_phases; count_i++)
//         {
//             GBs[cell]+=PopBal[active_phases[count_i]][cell]*PopBal[active_phases[count_i]][cell];

//             scalar initial_sum = Foam::pow(PopBal[active_phases[count_i]][cell],3.0)-PopBal[active_phases[count_i]][cell]+(2.0*PopBal[active_phases[count_i]][cell]*Foam::pow(1.0-xi[cell],2.0));

//                 Sum_Phases[cell]+=PopBal[active_phases[count_i]][cell];
//         for (label count_j = 0; count_j < n_active_phases; count_j++)
//         {
//                 if(active_phases[count_i] != active_phases[count_j]){
//                         initial_sum+=2.0*gamma_pf.value()*PopBal[active_phases[count_i]][cell]*Foam::pow(PopBal[active_phases[count_j]][cell],2.0);
//                 }
//         }

//         Sum_parts[active_phases[count_i]][cell]=initial_sum;


//         }



// }



forAll(mesh.C(), cell)
{
    const scalar xiCell = xi[cell];

    // --- First pass O(N): find active phases, accumulate totals ---
    scalar sumNiSq = 0.0;

    for (label count = 0; count < N_Ori; count++)
    {
        if (PopBal_av[count][cell] > SMALL)
        {
            const scalar ni = PopBal[count][cell];
            sumNiSq += ni*ni;
            Sum_Phases[cell] += ni;
        }
    }

    GBs[cell] = sumNiSq;  // GBs is just the sum of ni^2

    // --- Second pass O(N): compute Sum_parts using the identity ---
    //     cross term for phase i: Σ_{j≠i} 2γ·ni·nj²
    //                           = 2γ·ni·(sumNiSq - ni²)
    for (label count = 0; count < N_Ori; count++)
    {
        if (PopBal_av[count][cell] > SMALL)
        {
            const scalar ni = PopBal[count][cell];
            const scalar oneMinusXi = 1.0 - xiCell;

            Sum_parts[count][cell] =
                ni*ni*ni - ni
              + 2.0*ni*oneMinusXi*oneMinusXi
              + 2.0*gamma_pf.value()*ni*(sumNiSq - ni*ni);
        }
    }
}




Info << "Before solveNi " << endl;

for(label count = 0; count < N_Ori; count++)
{
    if (niActive[count])
    {
        volScalarField& ni = PopBal[count];
        volScalarField& Su=Sum_parts[count];
        
        solve
        (
            fvm::ddt(ni) 
        ==
            -Lg*m_g*Su
            +Lg*fvm::laplacian(Kg,ni)
        );
    }
}

phi_PF=0.5*(1.0-(Foam::tanh(1e5*((T/TLiquidus)-1.0))));

Info << "Before solve xi " << endl;

solve
(
    fvm::ddt(xi)
==
    -Lp*m_p*((2.0*xi*(1.0-phi_PF))-(2.0*phi_PF*(1.0-xi))) + Lp*m_g*(2.0*(1.0-xi)*GBs)
    +Lp*fvm::laplacian(Kp,xi)
);

xi=max(min(xi,1.0),0.0);

}
