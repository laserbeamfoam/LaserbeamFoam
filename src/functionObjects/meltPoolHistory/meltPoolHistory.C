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

    meltHistoryPtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "meltHistory",
                time_.timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
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
                time_.timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
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
                time_.timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
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
    previousLiquidFractionPtr_(),
    meltHistoryPtr_(),
    solidificationTimePtr_(),
    gradTSolPtr_()
{
    read(dict);
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
