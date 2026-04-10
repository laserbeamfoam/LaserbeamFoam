/*---------------------------------------------------------------------------*\
License
    This file is part of laserbeamFoam.

    laserbeamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    laserbeamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with laserbeamFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::mthdModel

Description
    Optional magneto-thermo-hydrodynamic model for the laserbeamFoam solvers.

\*---------------------------------------------------------------------------*/

#ifndef mthdModel_H
#define mthdModel_H

#include "pimpleControl.H"
#include "surfaceFields.H"
#include "volFields.H"

namespace Foam
{

class mthdModel
{
    // Private data

        fvMesh& mesh_;
        pimpleControl bpiso_;
        surfaceScalarField phiH_;
        volVectorField H_;
        volScalarField pH_;
        volVectorField B_;
        volVectorField J_;
        volScalarField divB_;
        volScalarField eleccond_;
        volScalarField magperm_;
        volScalarField oneoversigma_;
        volScalarField Joulesource_;
        volVectorField JcrossB_;


    // Private Member Functions

        static dimensionedScalar lookupElectricalConductivity
        (
            const dictionary& targetDictionary
        );

public:

    TypeName("mthdModel");

    // Constructors

        mthdModel(fvMesh& mesh, const word& propertiesDictName);

    // Destructor

        virtual ~mthdModel()
        {}


    // Static Member Functions

        static Switch electromagneticsEnabled
        (
            const fvMesh& mesh,
            const word& propertiesDictName
        );

        static dimensionedScalar readElectricalConductivity
        (
            const dictionary& targetDictionary
        );

        const volScalarField& Joulesource() const
        {
            return Joulesource_;
        }

        const volVectorField& JcrossB() const
        {
            return JcrossB_;
        }

        void correctIncompressibleProperties
        (
            const volScalarField& alpha1,
            const dictionary& transportPropertiesMetal,
            const dictionary& transportPropertiesGas
        );

        void resetTransportProperties();

        void accumulateTransportProperties
        (
            const volScalarField& alpha,
            const dimensionedScalar& electricalConductivity,
            const dimensionedScalar& magneticPermeability
        );

        void finaliseTransportProperties();

        void updateLorentzForce();

        void solve
        (
            const surfaceScalarField& phi,
            const volVectorField& U
        );
};

} // End namespace Foam

#endif
