if (!laser.radialPolarHeatSource())
{
    const faceList& ff = mesh.faces();
    const pointField& pp = mesh.points();
    const vectorField& CI = mesh.C();
    scalarField& yDimI = laser.yDim();

    forAll(CI, celli)
    {
        const cell& cc = mesh.cells()[celli];
        labelList pLabels(cc.labels(ff));
        pointField pLocal(pLabels.size(), vector::zero);

        forAll(pLabels, pointi)
        {
            pLocal[pointi] = pp[pLabels[pointi]];
        }

        yDimI[celli] =
            Foam::max(pLocal & vector(0, 1, 0))
          - Foam::min(pLocal & vector(0, 1, 0));
    }

    laser.yDim().correctBoundaryConditions();
}
