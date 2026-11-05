/*---------------------------------------------------------------------------*\
License
    This file is part of laserbeamFoam.

    laserbeamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    laserbeamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with laserbeamFoam.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::viscosityModels::temperatureTableViscosity

Description
    Kinematic viscosity nu(T) for the incompressible mixture, selected with
    transportModel temperatureTable; nu is any Function1 of T in the
    temperatureTableCoeffs sub-dictionary and is re-evaluated from the
    temperature field in correct().

\*---------------------------------------------------------------------------*/

#ifndef temperatureTableViscosity_H
#define temperatureTableViscosity_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "temperatureTable.H"

namespace Foam
{
namespace viscosityModels
{

class temperatureTableViscosity
:
    public viscosityModel
{
    // Private data

        dictionary coeffs_;

        // name of the temperture field
        word TName_;

        // used if the temperature field isnt registered yet
        scalar Tinit_;

        Foam::temperatureTable nuTable_;

        volScalarField nu_;


    // Private Member Functions

        void calcNu();


public:

    //- Runtime type information
    TypeName("temperatureTable");


    // Constructors

        temperatureTableViscosity
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~temperatureTableViscosity() = default;


    // Member Functions

        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        virtual void correct()
        {
            calcNu();
        }

        virtual bool read(const dictionary& viscosityProperties);

        const Foam::temperatureTable& nuTable() const
        {
            return nuTable_;
        }
};

} // End namespace viscosityModels
} // End namespace Foam

#endif

// ************************************************************************* //
