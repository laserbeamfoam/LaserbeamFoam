/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
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

\*---------------------------------------------------------------------------*/

#include "meltPoolHistory.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(meltPoolHistory, 0);
    addToRunTimeSelectionTable(functionObject, meltPoolHistory, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::meltPoolHistory::initialiseFields()
{
    if (previousLiquidFractionPtr_.valid())
    {
        return;
    }

    const volScalarField& liquidFraction =
        mesh_.lookupObject<volScalarField>(liquidFractionName_);

    // The fields are read from the time the run started from, so that a run
    // that continues from a written time goes on with its history instead of
    // starting it again at zero. They are written to the current time (the
    // instance is updated when they are written).
    const word& startTimeName = startTimeName_;

    previousLiquidFractionPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "previousLiquidFraction",
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            liquidFraction
        )
    );

    // Written with every time the solver writes (also with writeNow), not only
    // at the write times of this function object, so that a run can continue
    // from any written time
    meltHistoryPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "meltHistory",
                startTimeName,
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("meltHistory", dimless, 0.0)
        )
    );

    solidificationTimePtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "solidificationTime",
                startTimeName,
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("minusOne", dimTime, -1.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    gradTSolPtr_.reset
    (
        new volVectorField
        (
            IOobject
            (
                "gradTSol",
                startTimeName,
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimTemperature/dimLength, vector::zero),
            zeroGradientFvPatchVectorField::typeName
        )
    );
}


void Foam::functionObjects::meltPoolHistory::update()
{
    initialiseFields();

    const volScalarField& T = mesh_.lookupObject<volScalarField>
    (
        temperatureName_
    );
    const volScalarField& liquidFraction =
        mesh_.lookupObject<volScalarField>(liquidFractionName_);
    const volScalarField& alphaMetal =
        mesh_.lookupObject<volScalarField>(alphaMetalName_);

    volScalarField& previousLiquidFraction = previousLiquidFractionPtr_();
    volScalarField& meltHistory = meltHistoryPtr_();
    volScalarField& solidificationTime = solidificationTimePtr_();
    volVectorField& gradTSol = gradTSolPtr_();

    tmp<volVectorField> tgradT(fvc::grad(T));
    const volVectorField& gradT = tgradT();

    forAll(liquidFraction, cellI)
    {
        const scalar lf = liquidFraction[cellI];
        const scalar oldLf = previousLiquidFraction[cellI];

        if (alphaMetal[cellI] >= metalThreshold_ && lf >= liquidMetalThreshold_)
        {
            meltHistory[cellI] += 1.0;
        }

        if
        (
            lf <= solidificationThreshold_
         && oldLf > solidificationThreshold_
         && alphaMetal[cellI] > metalThreshold_ - SMALL
        )
        {
            scalar tSol = time_.deltaT().value();

            if (mag(oldLf - lf) > SMALL)
            {
                tSol = time_.deltaT().value()*oldLf/(oldLf - lf);
            }

            if (tSol < 0.0)
            {
                tSol = 0.0;
            }
            else if (tSol > time_.deltaT().value())
            {
                tSol = time_.deltaT().value();
            }

            solidificationTime[cellI] =
                time_.value() - time_.deltaT().value() + tSol;
            gradTSol[cellI] = gradT[cellI];
        }

        if (lf > solidificationThreshold_)
        {
            solidificationTime[cellI] = -1.0;
            gradTSol[cellI] = vector::zero;
        }
    }

    previousLiquidFraction = liquidFraction;

    meltHistory.correctBoundaryConditions();
    solidificationTime.correctBoundaryConditions();
    gradTSol.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::meltPoolHistory::meltPoolHistory
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    temperatureName_("T"),
    liquidFractionName_("epsilon1"),
    alphaMetalName_("alpha.metal"),
    metalThreshold_(0.5),
    liquidMetalThreshold_(0.5),
    solidificationThreshold_(SMALL),
    startTimeName_(runTime.timeName()),
    previousLiquidFractionPtr_(),
    meltHistoryPtr_(),
    solidificationTimePtr_(),
    gradTSolPtr_()
{
    read(dict);

    // Created now, at the time the run starts from, when the solver fields
    // exist already, so that a run that stops after its first time step still
    // writes them
    if (mesh_.foundObject<volScalarField>(liquidFractionName_))
    {
        initialiseFields();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::functionObjects::meltPoolHistory::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    temperatureName_ = dict.getOrDefault<word>("temperature", "T");
    liquidFractionName_ =
        dict.getOrDefault<word>("liquidFraction", "epsilon1");
    alphaMetalName_ = dict.getOrDefault<word>("alphaMetal", "alpha.metal");

    metalThreshold_ = dict.getOrDefault<scalar>("metalThreshold", 0.5);
    liquidMetalThreshold_ =
        dict.getOrDefault<scalar>("liquidMetalThreshold", 0.5);
    solidificationThreshold_ =
        dict.getOrDefault<scalar>("solidificationThreshold", SMALL);

    return true;
}


bool Foam::functionObjects::meltPoolHistory::execute()
{
    update();

    return true;
}


bool Foam::functionObjects::meltPoolHistory::write()
{
    initialiseFields();

    Log << type() << " " << name() << " write:" << nl
        << "    writing meltHistory" << nl
        << "    writing solidificationTime" << nl
        << "    writing gradTSol" << endl;

    meltHistoryPtr_().write();
    solidificationTimePtr_().write();
    gradTSolPtr_().write();

    return true;
}


// ************************************************************************* //
