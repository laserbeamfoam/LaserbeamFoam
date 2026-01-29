/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "compactRay.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(compactRay, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


compactRay::compactRay()
:
    position_(vector::zero),
    direction_(vector::one),
    power_(0),
    currentCell_(-1), 
    bounceCount_(0),
    globalRayIndex_(-1),
    active_(false),
    path_()
{}


compactRay::compactRay(const point& position, const vector& dir, scalar power)
:
    position_(position),
    direction_(dir),
    power_(power), 
    currentCell_(-1),
    bounceCount_(0),
    globalRayIndex_(-1),
    active_(true),
    path_()
{
    direction_ /= mag(direction_) + VSMALL;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
