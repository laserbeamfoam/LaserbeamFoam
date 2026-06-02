/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version: 4.0
    \\  /    A nd           | Web: http://www.foam-extend.org
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    Initialise the alpha.metal field for a powder bed multiphase simulation
    based on a "locations" file containing particle positions and radii
    (e.g. exported from LIGGGHTS).

Authors
    Petar Cosic
    Philip Cardiff
    Gowthaman Parivendhan

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvCFD.H"
#include "IOobject.H"
#include "volFields.H"
#include "timeSelector.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "IOstreams.H"

#include <fstream>
#include <set>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word getFieldName
(
    const Foam::argList& args,
    const Foam::word& fieldName,
    const Foam::word& defaultName = "alpha.phase"
)
{

    if (args.found(fieldName))
    {
        Foam::ITstream istr(args.lookup(fieldName));
        Foam::word fName;
        istr >> fName;
        Info << "Phase " << defaultName << " name: " << fName << nl;
        return fName;
    }

    Info << "Phase " << defaultName << " name: " << defaultName << nl;
    return defaultName;

}


std::vector<vector> subdivideCell
(
    const fvMesh& mesh,
    const label cellID,
    const int divisions,
    const vector& cellSize
)
{
    std::vector<vector> subCenters;
    subCenters.reserve(divisions*divisions*divisions);

    const vectorField& C = mesh.C();
    const vector cellC = C[cellID];
    const vector subSize = cellSize / divisions;
    const vector minCorner = cellC - 0.5*cellSize;

    for (int i = 0; i < divisions; ++i)
    {
        for (int j = 0; j < divisions; ++j)
        {
            for (int k = 0; k < divisions; ++k)
            {
                vector subC =
                    minCorner
                  + vector
                    (
                        (i + 0.5)*subSize.x(),
                        (j + 0.5)*subSize.y(),
                        (k + 0.5)*subSize.z()
                    );
                subCenters.push_back(subC);
            }
        }
    }

    return subCenters;
}


vector getCellSize
(
    const argList& args,
    const fvMesh& mesh,
    label repCellI = 0
)
{
    vector cellSizeVec(0, 0, 0);

    if (args.found("cellSize"))
    {
        Foam::ITstream istr(args.lookup("cellSize"));
        istr >> cellSizeVec;
        Info << "Using cellSize from command line: " << cellSizeVec << nl;
        return cellSizeVec;
    }
    else
    {
        if (mesh.nCells() == 0)
        {
            FatalErrorInFunction
                << "Mesh has no cells!" << nl
                << exit(FatalError);
        }

        if (repCellI < 0 || repCellI >= mesh.nCells())
        {
            FatalErrorInFunction
                << "Invalid cell index repCellI = " << repCellI << nl
                << exit(FatalError);
        }

        const cell& c = mesh.cells()[repCellI];
        const vector& p0 = mesh.points()[c[0]];
        const vector& p1 = mesh.points()[c[1]];

        vector cellDim;
        cellDim.x() = fabs(p1.x() - p0.x());
        cellDim.y() = fabs(p1.y() - p0.y());
        cellDim.z() = fabs(p1.z() - p0.z());

        scalar cellSIZE = 0.0;

        if (cellDim.x() > VSMALL)
        {
            cellSIZE = cellDim.x();
        }
        else if (cellDim.y() > VSMALL)
        {
            cellSIZE = cellDim.y();
        }
        else if (cellDim.z() > VSMALL)
        {
            cellSIZE = cellDim.z();
        }
        else
        {
            FatalErrorInFunction
                << "Cell size is zero in all directions." << nl
                << exit(FatalError);
        }

        cellSizeVec = vector(cellSIZE, cellSIZE, cellSIZE);
        Info << "Representative cell size: " << cellSizeVec << nl;
        return cellSizeVec;
    }
}


void printProgressBar(label current, label total, int barWidth = 50)
{
    scalar progress = scalar(current) / total;
    int pos = int(barWidth*progress);

    std::cout << "\r[";
    for (int i = 0; i < barWidth; ++i)
    {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress*100.0) << " %" << std::flush;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption("subDivisions", "int", "set the integer subDivisions value");
    argList::addOption("cellSize", "vector", "set the subdivision cell size");
    argList::addBoolOption("compressible", "enable compressible mode");
    argList::addOption("alpha.phase1", "word", "name of phase 1 volScalarField");
    argList::addOption("alpha.phase2", "word", "name of phase 2 volScalarField");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    bool compressible = args.found("compressible");

    word alphaPhase1Name = getFieldName(args, "alpha.phase1", "alpha.metal");

    word alphaPhase2Name;

    if(compressible){
        alphaPhase2Name = getFieldName(args, "alpha.phase2", "alpha.air");
    }

    volScalarField alphaPhase1
    (
        IOobject
        (
            alphaPhase1Name,
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0),
        "zeroGradient"
    );

    volScalarField alphaPhase2
    (
        IOobject
        (
            alphaPhase2Name,
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0),
        "zeroGradient"
    );

    label subDivisions = args.getOrDefault("subDivisions", 10);
    Info << "Subdivision subDivisions = " << subDivisions << nl;

    vector cellSize = getCellSize(args, mesh, 0);
    Info << "Using cell size for subdivision: " << cellSize << nl;

    // Read particle data
    Info << nl << "Reading the particle coordinates" << endl;

    std::ifstream loc(runTime.path() / "constant" / "location");
    if (!loc)
    {
        FatalError << nl << "Unable to open file 'location'" << nl
                   << exit(FatalError);
    }

    label count = 0, number = 0;
    scalarField particleR, particleX, particleY, particleZ;
    std::string line;

    while (std::getline(loc, line))
    {
        if (count == 3)
        {
            std::stringstream stream(line);
            stream >> number;

            particleR.setSize(number, 0.0);
            particleX.setSize(number, 0.0);
            particleY.setSize(number, 0.0);
            particleZ.setSize(number, 0.0);
        }
        else if (count > 8)
        {
            std::stringstream stream(line);
            stream >> particleX[count - 9]
                   >> particleY[count - 9]
                   >> particleZ[count - 9]
                   >> particleR[count - 9];
        }
        count++;
    }
    loc.close();

    // Read bedPlateDict
    const IOdictionary bedPlateDict
    (
        IOobject
        (
            "bedPlateDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const bool bedStatus(bedPlateDict.lookupOrDefault<bool>("Bed", false));
    scalar zmax = 0;

    if (bedStatus)
    {
        Info << "Reading Bed Plate properties" << endl;
        zmax = readScalar(bedPlateDict.lookup("zmax"));
    }
    else
    {
        Info << "Bed plate is switched off" << endl;
    }

    Info << "Setting field region values" << nl
         << "Number of particles in domain = " << number << endl;

    // Determine maximum height
    scalar maxHeight = -GREAT;
    label maxParticleIndex = -1;

    forAll(particleR, i)
    {
        scalar particleTopHeight = particleZ[i] + particleR[i];
        if (particleTopHeight > maxHeight)
        {
            maxHeight = particleTopHeight;
            maxParticleIndex = i;
        }
    }

    if (maxParticleIndex >= 0)
    {
        Info << "Highest particle at z = " << particleZ[maxParticleIndex]
             << " radius = " << particleR[maxParticleIndex]
             << " (total height: " << maxHeight << ")" << endl;
    }

    scalar layerHeight = maxHeight*1.0001;
    const vectorField& CI = mesh.C();
    scalarField& alphaPhase1I = alphaPhase1;

    forAll(CI, cellI)
    {
        if (cellI % ((CI.size() / 100 == 0) ? 1 : (CI.size() / 100)) == 0)
        {
            printProgressBar(cellI, CI.size());
        }

        const vector& curC = CI[cellI];
        scalar cellZ = curC.component(vector::Z);

        if (bedStatus && cellZ < zmax)
        {
            alphaPhase1I[cellI] = 1.0;
            continue;
        }

        if (cellZ >= zmax && cellZ <= layerHeight)
        {
            std::vector<vector> subCellCenters =
                subdivideCell(mesh, cellI, subDivisions, cellSize);

            for (Foam::label subCellI = 0; subCellI < Foam::label(subCellCenters.size()); ++subCellI) {

                const vector& subC = subCellCenters[subCellI];

                forAll(particleR, particleI)
                {
                    const scalar distance =
                        Foam::sqrt
                        (
                            Foam::pow(subC.x() - particleX[particleI], 2)
                          + Foam::pow(subC.y() - particleY[particleI], 2)
                          + Foam::pow(subC.z() - particleZ[particleI], 2)
                        );

                    if (distance <= particleR[particleI])
                    {
                        alphaPhase1I[cellI] += 1.0 / subCellCenters.size();
                        break;
                    }
                }
            }
        }
    }

    alphaPhase1.correctBoundaryConditions();
    Info << "\nWriting " << alphaPhase1.name()
         << " to time = " << runTime.timeName() << endl;
    alphaPhase1.write();

    if (compressible)
    {
        alphaPhase2 = 1 - alphaPhase1;
        alphaPhase2.write();
        Info << "\nWriting " << alphaPhase2.name()
         << " to time = " << runTime.timeName() << endl;
    }

    Info << "\nEnd" << endl;
    return 0;
}

// ************************************************************************* //
