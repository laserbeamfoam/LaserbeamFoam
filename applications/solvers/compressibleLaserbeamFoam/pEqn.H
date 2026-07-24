{

     if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }


    // volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rho*rAU())*fvc::ddtCorr(U, phi)
    );

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    PtrList<fvScalarMatrix> p_rghEqnComps(mixture.phases().size());

    label phasei = 0;
    forAllConstIters(mixture.phases(), phase)
    {
        const rhoThermo& thermo = phase().thermo();
        const tmp<volScalarField> trho(thermo.rho());
        const volScalarField& rho = trho();

        p_rghEqnComps.set
        (
            phasei,
            (
                fvc::ddt(rho) + thermo.psi()*correction(fvm::ddt(p_rgh))
              + fvc::div(phi, rho) - fvc::Sp(fvc::div(phi), rho)
            ).ptr()
        );

        ++phasei;
    }

    // Cache p_rgh prior to solve for density update
    volScalarField p_rgh_0(p_rgh);

    while (pimple.correctNonOrthogonal())
    {
        // Phase change, backward-Euler in pressure.  vDot was evaluated
        // (and clamped) once per step at start-of-step p; vDotP holds
        // -d(vDot)/dp >= 0 in the cells where no clamp bit.  The pair
        //     + fvm::Sp(vDotP, p_rgh) - vDotP*p_rgh.oldTime()
        // turns the frozen explicit source into
        //     vDot - vDotP*(p_rgh - p_rgh^n)
        //       ==  SUM_pairs sigma_e*A_i/delta*dV_i*(x_i*Psat_i - y_i*p_new)
        // in those cells: the solve itself relaxes each species' partial
        // pressure y_i*p toward its Raoult equilibrium x_i*Psat_i, the added
        // diagonal is positive (unconditionally stabilising), and it GROWS
        // as interface_thickness shrinks - inverting the explicit-lag
        // blow-up.  (interPhaseChangeFoam's vDotP pattern.)  For a pure
        // iteration stabiliser that leaves the converged source at the
        // frozen vDot, anchor on the current value instead: replace
        // p_rgh.oldTime() with p_rgh.
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
          - vDot
          + fvm::Sp(vDotP, p_rgh)
          - vDotP*p_rgh.oldTime()
        );

        tmp<fvScalarMatrix> p_rghEqnComp;

        phasei = 0;
        forAllConstIters(mixture.phases(), phase)
        {
            tmp<fvScalarMatrix> hmm
            (
                (max(phase(), scalar(0))/phase().thermo().rho())
               *p_rghEqnComps[phasei]
            );

            if (phasei == 0)
            {
                p_rghEqnComp = hmm;
            }
            else
            {
                p_rghEqnComp.ref() += hmm;
            }

            ++phasei;
        }

        solve
        (
            p_rghEqnComp
          + p_rghEqnIncomp,
            p_rgh.select(pimple.finalInnerIter())
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phasei = 0;
            for (phaseModel& phase : mixture.phases())
            {
                phase.dgdt() =
                    pos0(phase)
                  *(p_rghEqnComps[phasei] & p_rgh)/phase.thermo().rho();

                ++phasei;
            }

            phi = phiHbyA + p_rghEqnIncomp.flux();

            U = HbyA
              + rAU()*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);
            U.correctBoundaryConditions();
        }
    }

    p = max(p_rgh + mixture.rho()*gh, pMin);

    // Update densities from change in p_rgh
    mixture.correctRho(p_rgh - p_rgh_0);
    rho = mixture.rho();

    // Correct p_rgh for consistency with p and the updated densities
    p_rgh = p - rho*gh;
    p_rgh.correctBoundaryConditions();

    K = 0.5*magSqr(U);

    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min(p_rgh) " << min(p_rgh).value() << endl;
    Info<< "min(T) " << min(T).value() << endl;

    if (!correctPhi)
    {
        rAU.clear();
    }

}