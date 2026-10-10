/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    TesselateFoam

Description
    Utility to set phase fields using Voronoi tessellation.
    Each mesh cell is assigned to the nearest seed point (nearest-neighbour
    assignment is equivalent to Voronoi tessellation).

Author
    Tom Flint

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Read xi field
    volScalarField xi
    (
        IOobject
        (
            "xi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("xi", dimensionSet(0, 0, 0, 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    // Read phase field properties dictionary
    IOdictionary PhaseFieldProperties
    (
        IOobject
        (
            "PhaseFieldProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    label N_Seeds(readLabel(PhaseFieldProperties.lookup("N_Seeds")));
    label N_Ori = N_Seeds + 1;

    Info << "Number of SEEDS: " << N_Seeds << endl;

    // Create phase field order parameters
    PtrList<volScalarField> PopBal(N_Ori);

    for (label count = 0; count < N_Ori; count++)
    {
        word name_ni("n." + name(count));

        PopBal.set
        (
            count,
            new volScalarField
            (
                IOobject
                (
                    name_ni,
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("ni", dimensionSet(0, 0, 0, 0, 0), 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }

    // Read tessellation bounds
    scalar Xmin(readScalar(PhaseFieldProperties.lookup("Xmin")));
    scalar Xmax(readScalar(PhaseFieldProperties.lookup("Xmax")));
    scalar Ymin(readScalar(PhaseFieldProperties.lookup("Ymin")));
    scalar Ymax(readScalar(PhaseFieldProperties.lookup("Ymax")));
    scalar Zmin(readScalar(PhaseFieldProperties.lookup("Zmin")));
    scalar Zmax(readScalar(PhaseFieldProperties.lookup("Zmax")));

    // -----------------------------------------------------------------------
    // Generate random seed points within the bounding box
    // -----------------------------------------------------------------------
    Random rndGen(1234); // Fixed seed for reproducibility

    List<vector> seedPoints(N_Seeds);

    for (label i = 0; i < N_Seeds; i++)
    {
        scalar x = Xmin + rndGen.sample01<scalar>() * (Xmax - Xmin);
        scalar y = Ymin + rndGen.sample01<scalar>() * (Ymax - Ymin);
        scalar z = Zmin + rndGen.sample01<scalar>() * (Zmax - Zmin);
        seedPoints[i] = vector(x, y, z);

        Info << "Seed[" << i << "] = " << seedPoints[i] << endl;
    }

    // -----------------------------------------------------------------------
    // Voronoi tessellation by nearest-neighbour assignment
    // For each mesh cell, find the closest seed and set that field to 1.0
    // -----------------------------------------------------------------------
    Info << "Assigning cells to nearest seed point..." << endl;

    forAll(mesh.C(), celli)
    {
        // Only assign cells where the last PopBal field is not already set
        // and xi indicates liquid (matching original logic)
        if
        (
            PopBal[N_Seeds][celli] < 1.0 - SMALL
         && xi[celli] > 1.0 - SMALL
        )
        {
            const vector& cellCentre = mesh.C()[celli];

            // Find nearest seed
            label nearestSeed = -1;
            scalar minDistSq = GREAT;

            forAll(seedPoints, si)
            {
                scalar distSq = magSqr(cellCentre - seedPoints[si]);

                if (distSq < minDistSq)
                {
                    minDistSq = distSq;
                    nearestSeed = si;
                }
            }

            if (nearestSeed >= 0)
            {
                PopBal[nearestSeed][celli] = 1.0;
            }
        }
    }

    // Write all phase fields
    for (label count = 0; count < N_Ori; count++)
    {
        PopBal[count].write();
    }

    // -----------------------------------------------------------------------
    // Populate grainNum field
    // -----------------------------------------------------------------------
    volScalarField grainNum
    (
        IOobject
        (
            "grainNum",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("grainNum", dimensionSet(0, 0, 0, 0, 0), -1.0),
        zeroGradientFvPatchScalarField::typeName
    );

    volScalarField maxNiVal
    (
        IOobject
        (
            "maxNiVal",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("maxNiVal", dimensionSet(0, 0, 0, 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    scalar grainNumThreshold = 0.5;

    forAll(grainNum, i)
    {
        forAll(PopBal, j)
        {
            if
            (
                (PopBal[j][i] > grainNumThreshold)
             && (PopBal[j][i] > maxNiVal[i])
            )
            {
                maxNiVal[i] = PopBal[j][i];
                grainNum[i] = j;
            }
        }
    }

    grainNum.write();

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
