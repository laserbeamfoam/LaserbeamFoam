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
    // ddtCorr weight and gravity split use rhoM to match UEqn.H; compressibility below stays EOS
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rhoM*rAU())*fvc::ddtCorr(U, phi)
    );

    surfaceScalarField phig
    (
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rhoM)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf);

    // per-phase (1/rho_i)*D(rho_i)/Dt: only feeds phase.dgdt() for the beta transport now
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

    // adiabatic compressibility at fixed composition, chi*sum_i alpha_i*psi_i/rho_i [1/Pa]
    volScalarField psiBar
    (
        IOobject("psiBar", runTime.timeName(), mesh),
        mesh,
        dimensionedScalar(dimless/dimPressure, Zero)
    );
    for (const phaseModel& phase : mixture.phases())
    {
        psiBar += max(phase, scalar(0))*phase.thermo().psi()/phase.thermo().rho();
    }
    psiBar *= chiClosure;

    const dimensionedScalar rDeltaTp(1.0/runTime.deltaT());

    // Cache p_rgh prior to solve for density update
    volScalarField p_rgh_0(p_rgh);

    while (pimple.correctNonOrthogonal())
    {
        // closure source (heating, phase change, transport residual) plus the
        // backward-Euler pressure coupling of the phase-change rate
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
          - closureSrc
          + fvm::Sp(vDotP, p_rgh)
          - vDotP*p_rgh.oldTime()
        );

        // the volume constraint linearised about the closure state:
        // chi*psi*(p - p_S)/dt + div(phi) = chi*(s - 1)/dt   (Zenz: Psi*(p - p*))
        fvScalarMatrix p_rghEqnComp
        (
            fvm::Sp(psiBar*rDeltaTp, p_rgh)
          - psiBar*rDeltaTp*p_rghS
        );

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

    // hydrostatic split with rhoM (frozen over the correctors, so the round trip is exact)
    p = max(p_rgh + rhoM*gh, pMin);

    // Update densities from change in p_rgh
    mixture.correctRho(p_rgh - p_rgh_0);
    rho = mixture.rho();

    // Correct p_rgh for consistency with p
    p_rgh = p - rhoM*gh;
    p_rgh.correctBoundaryConditions();

    K = 0.5*magSqr(U);

    fvc::correctUf(Uf, U, phi);

    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min,max(p_rgh) " << min(p_rgh).value()
        << " " << max(p_rgh).value() << endl;
    Info<< "min(T) " << min(T).value() << endl;
    Info<< "max|vDot|*dt (phase change) " << max(mag(vDot)).value()*mesh.time().deltaTValue()
        << ", max|closureSrc|*dt " << max(mag(closureSrc)).value()*mesh.time().deltaTValue()
        << endl;

    if (!correctPhi)
    {
        rAU.clear();
    }

}
