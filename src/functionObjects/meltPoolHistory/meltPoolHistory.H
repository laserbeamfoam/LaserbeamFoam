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

    The fields meltHistory, solidificationTime and gradTSol are written at
    every write time of the solver (including writeNow), whatever the
    writeControl of the function object, and read from the start time, so
    that a run continued from any written time keeps its history.

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

        //- Name of the temperature field
        word temperatureName_;

        //- Name of the liquid fraction field
        word liquidFractionName_;

        //- Name of the metal volume fraction field
        word alphaMetalName_;

        //- Metal volume fraction above which a cell is metal
        scalar metalThreshold_;

        //- Liquid fraction above which a metal cell is liquid
        scalar liquidMetalThreshold_;

        //- Liquid fraction below which a cell is solid
        scalar solidificationThreshold_;

        //- Name of the time the run started from, where the fields are read
        word startTimeName_;

        //- Liquid fraction at the previous time step
        autoPtr<volScalarField> previousLiquidFractionPtr_;

        //- Number of time steps each metal cell has been liquid
        autoPtr<volScalarField> meltHistoryPtr_;

        //- Time of the last solidification of each cell (-1 if none)
        autoPtr<volScalarField> solidificationTimePtr_;

        //- Temperature gradient at the last solidification of each cell
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
