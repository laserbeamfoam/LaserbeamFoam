/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "phaseModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::phaseModel::phaseState Foam::phaseModel::readPhaseState
(
    const dictionary& dict
) const
{
    word stateWord = dict.getOrDefault<word>("phaseState", "condensed");
    
    if (stateWord == "condensed")
    {
        return phaseState::CONDENSED;
    }
    else if (stateWord == "gaseous" || stateWord == "gas" || stateWord == "vapour")
    {
        return phaseState::GASEOUS;
    }
    else if (stateWord == "solid")
    {
        return phaseState::SOLID;
    }
    else if (stateWord == "plasma")
    {
        return phaseState::PLASMA;
    }
    else
    {
        WarningInFunction
            << "Unknown phase state '" << stateWord 
            << "' for phase " << name_
            << ". Valid options: condensed, gaseous, solid, plasma"
            << nl << "Defaulting to condensed."
            << endl;
        
        return phaseState::CONDENSED;
    }
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModel::phaseModel
(
    const word& phaseName,
    const volScalarField& p,
    const volScalarField& T
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            p.mesh().time().timeName(),
            p.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        p.mesh()
    ),
    name_(phaseName),
    p_(p),
    T_(T),
    thermo_(nullptr),
    dgdt_
    (
        IOobject
        (
            IOobject::groupName("dgdt", phaseName),
            p.mesh().time().timeName(),
            p.mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        p.mesh(),
        dimensionedScalar(dimless/dimTime, Zero)
    ),
    state_(phaseState::CONDENSED)
{
    {
        volScalarField Tp(IOobject::groupName("T", phaseName), T);
        Tp.write();
    }

    thermo_ = rhoThermo::New(p.mesh(), phaseName);
    thermo_->validate(phaseName, "e");

        // Read phase state from thermophysicalProperties dictionary
    IOdictionary phaseDictionary
    (
        IOobject
        (
            "thermophysicalProperties." + phaseName,
            p.mesh().time().constant(),
            p.mesh(),
            IOobject::MUST_READ_IF_MODIFIED
        )
    );
    
    state_ = readPhaseState(phaseDictionary);
    
    Info<< "Phase " << phaseName << " state: " << stateAsWord() << endl;




    correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseModel> Foam::phaseModel::clone() const
{
    NotImplemented;
    return nullptr;
}


void Foam::phaseModel::correct()
{
    thermo_->he() = thermo_->he(p_, T_);
    thermo_->correct();
}

Foam::word Foam::phaseModel::stateAsWord() const
{
    switch (state_)
    {
        case phaseState::CONDENSED:
            return "condensed";
        case phaseState::GASEOUS:
            return "gaseous";
        case phaseState::SOLID:
            return "solid";
        case phaseState::PLASMA:
            return "plasma";
        default:
            return "unknown";
    }
}

// ************************************************************************* //
