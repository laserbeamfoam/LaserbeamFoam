{

     if (correctPhi)
    {
        rAU.ref() = 1.0/UEqn.A();
    }
    else
    {
        rAU = 1.0/UEqn.A();
    }


    // ddtCorr weight and gravity split use rhoM to match UEqn.H; compressibility below stays EOS
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU()));
    volVectorField HbyA(constrainHbyA(rAU()*UEqn.H(), U, p_rgh));
    // ddtCorr of Choi 1999 (https://doi.org/10.1080/104077899274679) written out: OpenFOAM's empirical
    // weight 1 - |phi0 - flux(U0)|/|phi0| (https://www.tfd.chalmers.se/~hani/kurser/OS_CFD_2022/lectureNotes/ddtCorr.pdf)
    // reads a phase-change flux jump as an error (odd-even p, test 09): phase-change faces keep weight 1
    tmp<surfaceScalarField> tphi0
    (
        mesh.dynamic()
      ? (mesh.Sf() & Uf().oldTime())        // as ddtCorr(U, phi, Uf): AMR maps phi_0 face by face (test 12)
      : tmp<surfaceScalarField>(phi.oldTime())
    );
    const surfaceScalarField& phi0 = tphi0();
    // memory of the smooth part: the Stefan jump the old interpolation missed (the old step's
    // phiJump) is not memory, it is put back for this step below (stefanFlow.H)
    surfaceScalarField phiCorrDdt
    (
        "phiCorrDdt",
        phi0 - fvc::flux(U.oldTime()) - uJump.oldTime()*mesh.magSf()
    );
    surfaceScalarField ddtCoeff
    (
        "ddtCoeff",
        scalar(1)
      - min
        (
            mag(phiCorrDdt)
           /(mag(phi0) + dimensionedScalar("small", phi.dimensions(), SMALL)),
            scalar(1)
        )
    );
    {
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        const scalarField& vI = vDot.primitiveField();
        scalarField& cI = ddtCoeff.primitiveFieldRef();
        forAll(own, facei)
        {
            if (vI[own[facei]] != 0 || vI[nei[facei]] != 0) cI[facei] = 1.0;
        }

        surfaceScalarField::Boundary& cBf = ddtCoeff.boundaryFieldRef();
        forAll(U.boundaryField(), patchi)
        {
            const fvPatch& fvp = mesh.boundary()[patchi];
            if (U.boundaryField()[patchi].fixesValue() || isA<cyclicAMIFvPatch>(fvp))
            {
                cBf[patchi] = 0.0;
            }
            else if (fvp.coupled())
            {
                const scalarField vNb(vDot.boundaryField()[patchi].patchNeighbourField());
                const labelUList& fc = fvp.faceCells();
                forAll(fc, i)
                {
                    if (vI[fc[i]] != 0 || vNb[i] != 0) cBf[patchi][i] = 1.0;
                }
            }
        }
    }

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + phiJump                     // the Stefan jump, which the interpolation of HbyA misses
      + fvc::interpolate(rhoM*rAU())*ddtCoeff*phiCorrDdt/runTime.deltaT()
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
        // closure source: heating, phase change, transport residual
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
          - closureSrc
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
