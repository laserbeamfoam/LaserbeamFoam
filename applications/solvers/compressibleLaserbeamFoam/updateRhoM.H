// rhoM = sum_i c_i, floored at rhoMMin; once per time step, after mixture.solve()
{
    const tmp<volScalarField> trhoMledger(mixture.rhoM());
    const volScalarField& rhoMledger = trhoMledger();

    label nFloor = 0;
    {
        const scalarField& rI = rhoMledger.primitiveField();
        forAll(rI, celli)
        {
            if (rI[celli] < rhoMMin.value())
            {
                ++nFloor;
            }
        }
        reduce(nFloor, sumOp<label>());
    }

    rhoM = max(rhoMledger, rhoMMin);

    // EOS vs ledger density (1 everywhere = closure exact), on the true ledger
    {
        const scalarField& rhoI  = rho.primitiveField();
        const scalarField& rhoMI = rhoM.primitiveField();
        const scalarField& V     = mesh.V().field();

        const scalarField& CvI = CvMass.primitiveField();
        const scalarField& TI  = T.primitiveField();

        scalar ratioMin = GREAT;
        scalar ratioMax = -GREAT;
        label nOff = 0;
        scalar massEOS = 0.0;
        scalar massLedger = 0.0;
        scalar massGuard = 0.0;      // what ledgerFloorEOS would add
        scalar enthalpyGuard = 0.0;  // its sensible heat, untracked when the guard engages/releases

        forAll(rhoI, celli)
        {
            const scalar ratio = rhoI[celli]/rhoMI[celli];
            ratioMin = Foam::min(ratioMin, ratio);
            ratioMax = Foam::max(ratioMax, ratio);
            if (ratio < 0.5 || ratio > 2.0)
            {
                ++nOff;
            }
            massEOS    += rhoI[celli]*V[celli];
            massLedger += rhoMI[celli]*V[celli];
            const scalar dm = Foam::max(rhoI[celli] - rhoMI[celli], 0.0)*V[celli];
            massGuard     += dm;
            enthalpyGuard += dm*CvI[celli]*TI[celli];
        }
        reduce(ratioMin, minOp<scalar>());
        reduce(ratioMax, maxOp<scalar>());
        reduce(nOff, sumOp<label>());
        reduce(massEOS, sumOp<scalar>());
        reduce(massLedger, sumOp<scalar>());
        reduce(massGuard, sumOp<scalar>());
        reduce(enthalpyGuard, sumOp<scalar>());

        Info<< "rho(EOS)/rhoM(ledger): min = " << ratioMin
            << ", max = " << ratioMax
            << ", cells outside [0.5,2] = " << nOff
            << "; total mass EOS " << massEOS
            << " kg vs ledger " << massLedger << " kg"
            << "; guard mass " << massGuard << " kg ("
            << enthalpyGuard << " J)" << endl;
    }

    label nGuard = 0;
    if (ledgerFloorEOS)
    {
        // under-full cells carry the EOS mass the pressure equation assumes; count raises above 1%
        const scalarField& rhoI = rho.primitiveField();
        const scalarField& rhoMI = rhoM.primitiveField();
        forAll(rhoMI, celli)
        {
            if (rhoI[celli] > 1.01*rhoMI[celli])
            {
                ++nGuard;
            }
        }
        reduce(nGuard, sumOp<label>());
        rhoM = max(rhoM, rho);
    }
    rhoM.correctBoundaryConditions();

    if (nFloor || nGuard)
    {
        Info<< "rhoM floored at " << rhoMMin.value() << " kg/m3 in "
            << nFloor << " cell(s)";
        if (ledgerFloorEOS)
        {
            Info<< ", raised to EOS rho by >1% in " << nGuard << " cell(s)";
        }
        Info<< endl;
    }
}