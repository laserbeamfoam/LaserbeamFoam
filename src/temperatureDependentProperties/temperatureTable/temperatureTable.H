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
    Foam::temperatureTable

Description
    Temperature dependant material property: any Function1 of T (constant,
    table, tableFile, polynomial, ...) sampled once onto a uniform T grid so
    the per cell evalution is a single interpolation.

\*---------------------------------------------------------------------------*/

#ifndef temperatureTable_H
#define temperatureTable_H

#include "scalarField.H"
#include "string.H"
#include "dictionary.H"
#include "Function1.H"
#include "Polynomial.H"

namespace Foam
{

class temperatureTable
{
    // Private data

        word name_;

        string description_;

        scalar Tmin_;

        scalar Tmax_;

        scalar deltaT_;

        scalar rDeltaT_;

        // values at Tmin_ + i*deltaT_
        scalarField values_;


    // Private Member Functions

        void readControls(const dictionary& controls);

        label nSamples() const;

        void sample(const Function1<scalar>& f);

        void sample(const Polynomial<8>& p);


public:

    // Constructors

        // Function1 entry 'keyword' of dict, or the legacy Polynomial<8> entry
        // if the keyword is absent
        temperatureTable
        (
            const word& keyword,
            const dictionary& dict,
            const dictionary& controls,
            const word& legacyPolynomialKeyword = word::null
        );

        // constant property
        temperatureTable
        (
            const word& name,
            const scalar value,
            const dictionary& controls
        );


    // Member Functions

        const word& name() const { return name_; }

        const string& description() const { return description_; }

        scalar Tmin() const { return Tmin_; }

        scalar Tmax() const { return Tmax_; }

        scalar deltaT() const { return deltaT_; }

        label size() const { return values_.size(); }

        bool isConstant() const;

        // property at T, clamped to the cached range
        inline scalar value(const scalar T) const;

        void values(const UList<scalar>& T, UList<scalar>& result) const;

        tmp<scalarField> values(const UList<scalar>& T) const;

        // print the specifcation and the values at Tsample to Info
        void printInfo(const UList<scalar>& Tsample) const;
};

} // End namespace Foam

#include "temperatureTableI.H"

#endif

// ************************************************************************* //
