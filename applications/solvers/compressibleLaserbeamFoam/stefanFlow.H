// Stefan flow of this step's condensation (UEqn.H, pEqn.H). A smeared jump [u] = j[1/rho] gives the
// viscous term a spurious stress ~mu*[u]/dx (https://doi.org/10.1016/j.jcp.2024.113704), removed here
// as ghost-fluid methods remove it (https://doi.org/10.1016/j.jcp.2016.04.031). Evaporation: REVIEW.md
{
    const volScalarField& sGas = mixture.phaseChangeGasDilatation();
    const volScalarField& sCond = mixture.phaseChangeCondensedDilatation();

    // the donors: cells the library drew condensing vapour from this step
    volScalarField donor
    (
        IOobject
        (
            "condensationDonor",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    );
    forAll(donor, celli)
    {
        if (sGas[celli] < 0) donor[celli] = 1;
    }
    donor.correctBoundaryConditions();

    // dilatation of condensation: vapour volume drawn (< 0) plus condensate volume formed (> 0)
    const dimensionedScalar zeroRate(sGas.dimensions(), Zero);
    volScalarField vDotPC("vDotPC", min(sGas, zeroRate) + max(sCond, zeroRate));

    // closed domain: the net volume change is taken up by the expansion of the whole gas
    if (psiStefan.needReference())
    {
        const scalarField wGas(max(1.0 - condensate.primitiveField(), scalar(0)));
        const scalar netV = gSum(vDotPC.primitiveField()*mesh.V().field());
        const scalar gasV = gSum(wGas*mesh.V().field());

        if (gasV > VSMALL)
        {
            vDotPC.primitiveFieldRef() -= wGas*(netV/gasV);
        }
        else
        {
            vDotPC.primitiveFieldRef() -= netV/gSum(mesh.V().field());
        }
    }

    // potential flow through the gas, mobility 1 - condensate (floored, harmonic on faces): the
    // Stefan potential of Scapin et al. 2020, https://doi.org/10.1016/j.jcp.2020.109251
    const volScalarField kGas("kGas", max(1.0 - condensate, dimensionedScalar(dimless, 1e-3)));
    const surfaceScalarField kGasF("kGasF", 1.0/fvc::interpolate(1.0/kGas));

    fvScalarMatrix psiStefanEqn
    (
        fvm::laplacian(kGasF, psiStefan) == vDotPC
    );
    psiStefanEqn.setReference(Pstream::master() ? 0 : -1, 0.0);

    // p_rgh solver settings unless fvSolution has a psiStefan entry
    const dictionary& fvSolSolvers =
        static_cast<const fvSolution&>(mesh).subDict("solvers");

    psiStefanEqn.solve
    (
        fvSolSolvers.found(psiStefan.name())
      ? fvSolSolvers.subDict(psiStefan.name())
      : mesh.solverDict(p_rgh.select(true))
    );

    const surfaceScalarField phiStefan("phiStefan", psiStefanEqn.flux());
    Upc = fvc::reconstruct(phiStefan);

    // condensing vapour leaves a donor with the donor's velocity (library momentum pairing), so the
    // donor keeps the arrival velocity and the whole jump sits on its face to the condensed side:
    // Upc += s*h/2 towards the gas, h the cell's thickness along the interface normal
    {
        const volVectorField gradBeta(fvc::grad(condensate));
        const vectorField& Sa = mesh.faceAreas();
        const cellList& cellFaces = mesh.cells();
        const scalarField& Vc = mesh.V();
        vectorField& UpcI = Upc.primitiveFieldRef();

        forAll(UpcI, celli)
        {
            const scalar magGradBeta = mag(gradBeta[celli]);
            if (donor[celli] == 0 || magGradBeta < VSMALL) continue;

            const vector nGas(-gradBeta[celli]/magGradBeta);
            scalar aPerp = 0;
            for (const label facei : cellFaces[celli])
            {
                aPerp += 0.5*mag(Sa[facei] & nGas);
            }
            UpcI[celli] += 0.5*sGas[celli]*(Vc[celli]/aPerp)*nGas;
        }
    }
    Upc.correctBoundaryConditions();

    // U is mass-averaged and Upc is the gas volume flux per area, so the Stefan part of U is
    // (rho_gas/rho)*Upc: 1 in clean gas, ~0 where liquid carries the cell's momentum (an anchor)
    {
        scalarField cGas(mesh.nCells(), Zero);
        scalarField aGas(mesh.nCells(), Zero);
        scalarField cAll(mesh.nCells(), Zero);

        label phasei = 0;
        for (const phaseModel& ph : mixture.phases())
        {
            const scalarField cPh(max(mixture.c()[phasei].primitiveField(), scalar(0)));
            cAll += cPh;
            if (ph.isGaseous())
            {
                cGas += cPh;
                aGas += max(ph.primitiveField(), scalar(0));
            }
            ++phasei;
        }

        vectorField& UpcI = Upc.primitiveFieldRef();
        forAll(UpcI, celli)
        {
            UpcI[celli] *=
                (aGas[celli] > SMALL && cAll[celli] > VSMALL)
              ? min(cGas[celli]/aGas[celli]/cAll[celli], scalar(1))
              : scalar(0);
        }
    }
    Upc.correctBoundaryConditions();

    // the flux the interpolation of U's Stefan part misses on the donors' faces (zero elsewhere): the
    // face velocity jump of the ghost-fluid projection, https://doi.org/10.1016/j.jcp.2014.01.014
    phiJump = phiStefan - fvc::flux(Upc);
    {
        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();
        scalarField& jI = phiJump.primitiveFieldRef();

        forAll(jI, facei)
        {
            if (donor[own[facei]] == 0 && donor[nei[facei]] == 0) jI[facei] = 0;
        }

        forAll(phiJump.boundaryField(), patchi)
        {
            fvsPatchScalarField& jp = phiJump.boundaryFieldRef()[patchi];

            if (!jp.coupled())
            {
                jp = 0;
                continue;
            }

            const scalarField donorNb(donor.boundaryField()[patchi].patchNeighbourField());
            const labelUList& fc = mesh.boundary()[patchi].faceCells();
            forAll(jp, i)
            {
                if (donor[fc[i]] == 0 && donorNb[i] == 0) jp[i] = 0;
            }
        }
    }
    uJump = phiJump/mesh.magSf();

    // largest jump on any face, processor and cyclic faces included
    scalar maxJump = 0;
    for (const scalar u : uJump.primitiveField())
    {
        maxJump = max(maxJump, mag(u));
    }
    forAll(uJump.boundaryField(), patchi)
    {
        if (uJump.boundaryField()[patchi].coupled())
        {
            for (const scalar u : uJump.boundaryField()[patchi])
            {
                maxJump = max(maxJump, mag(u));
            }
        }
    }
    reduce(maxJump, maxOp<scalar>());

    Info<< "Stefan flow of condensation: max|Upc| " << gMax(mag(Upc.primitiveField())())
        << " m/s, max jump the interpolation misses " << maxJump << " m/s" << endl;
}
