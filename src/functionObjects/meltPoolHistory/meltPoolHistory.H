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

Class
    Foam::functionObjects::meltPoolHistory

Description
    Tracks melt-pool history fields for laserbeamFoam-style simulations.

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    meltPoolHistory
    {
        type            meltPoolHistory;
        libs            ("liblaserbeamFoamFunctionObjects.so");

        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;

        temperature     T;
        liquidFraction  epsilon1;
        alphaMetal      alpha.metal;

        liquidMetalThreshold      0.5;
        solidificationThreshold   1e-15;
    }
    \endverbatim

SourceFiles
    meltPoolHistory.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_meltPoolHistory_H
#define functionObjects_meltPoolHistory_H

#include "fvMeshFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                        Class meltPoolHistory Declaration
\*---------------------------------------------------------------------------*/

class meltPoolHistory
:
    public fvMeshFunctionObject
{
    // Private Data

        word temperatureName_;
        word liquidFractionName_;
        word alphaMetalName_;

        scalar metalThreshold_;
        scalar liquidMetalThreshold_;
        scalar solidificationThreshold_;

        autoPtr<volScalarField> previousLiquidFractionPtr_;
        autoPtr<volScalarField> meltHistoryPtr_;
        autoPtr<volScalarField> solidificationTimePtr_;
        autoPtr<volVectorField> gradTSolPtr_;


    // Private Member Functions

        void initialiseFields();

        void update();


public:

    //- Runtime type information
    TypeName("meltPoolHistory");


    // Constructors

        //- Construct from Time and dictionary
        meltPoolHistory
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        meltPoolHistory(const meltPoolHistory&) = delete;

        //- No copy assignment
        void operator=(const meltPoolHistory&) = delete;


    //- Destructor
    virtual ~meltPoolHistory() = default;


    // Member Functions

        //- Read the function object controls
        virtual bool read(const dictionary& dict);

        //- Update the melt-pool history fields
        virtual bool execute();

        //- Write the melt-pool history fields
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
