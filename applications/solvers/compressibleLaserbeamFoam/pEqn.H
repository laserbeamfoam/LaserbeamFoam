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
        // Phase change, backward-Euler in pressure.
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


            {
                const scalar UMax =
                    pimple.dict().getOrDefault<scalar>("UMax", 1000.0);

                vectorField& UI = U.primitiveFieldRef();
                label nUClip = 0;
                forAll(UI, celli)
                {
                    const scalar magU = mag(UI[celli]);
                    if (magU > UMax)
                    {
                        UI[celli] *= UMax/magU;
                        ++nUClip;
                    }
                }
                reduce(nUClip, sumOp<label>());

                if (nUClip)
                {
                    U.correctBoundaryConditions();
                    Info<< "pEqn: |U| ceiling (" << UMax
                        << " m/s) engaged in " << nUClip
                        << " cell(s)" << endl;
                }

                label nPhiClip = 0;
                const scalarField& magSfI = mesh.magSf().primitiveField();
                scalarField& phiI = phi.primitiveFieldRef();
                forAll(phiI, facei)
                {
                    const scalar phiMax = UMax*magSfI[facei];
                    if (mag(phiI[facei]) > phiMax)
                    {
                        phiI[facei] = sign(phiI[facei])*phiMax;
                        ++nPhiClip;
                    }
                }

                surfaceScalarField::Boundary& phiBf = phi.boundaryFieldRef();
                const surfaceScalarField::Boundary& magSfBf =
                    mesh.magSf().boundaryField();
                forAll(phiBf, patchi)
                {
                    scalarField& phiP = phiBf[patchi];
                    const scalarField& magSfP = magSfBf[patchi];
                    forAll(phiP, facei)
                    {
                        const scalar phiMax = UMax*magSfP[facei];
                        if (mag(phiP[facei]) > phiMax)
                        {
                            phiP[facei] = sign(phiP[facei])*phiMax;
                            ++nPhiClip;
                        }
                    }
                }
                reduce(nPhiClip, sumOp<label>());

                if (nPhiClip)
                {
                    Info<< "pEqn: |phi| ceiling engaged on " << nPhiClip
                        << " face(s)" << endl;
                }
            }
            // ----------------------------------------------------------------
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

    fvc::correctUf(Uf, U, phi);

    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min,max(p_rgh) " << min(p_rgh).value()
        << " " << max(p_rgh).value() << endl;
    Info<< "min(T) " << min(T).value() << endl;
    Info<< "max|vDot|*dt " << max(mag(vDot)).value()*mesh.time().deltaTValue()
        << endl;

    if (!correctPhi)
    {
        rAU.clear();
    }

}