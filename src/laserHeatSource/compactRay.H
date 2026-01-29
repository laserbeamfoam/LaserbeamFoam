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

Class
    compactRay

Description
    Class for storing information about individual rays propagating through
    a polMesh in a ray tracing algorithm.

Author
    Thomas Flint, The University of Manchester.
    Philip Cardiff, University College Dublin.

SourceFiles
    compactRay.C

\*---------------------------------------------------------------------------*/

#ifndef compactRay_H
#define compactRay_H

#include "point.H"
#include "vector.H"
#include "DynamicList.H"
#include "HashSet.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class compactRay Declaration
\*---------------------------------------------------------------------------*/


class compactRay
{
    public:

        // Public data

            //- Current position
            point position_;

            //- Current direction
            vector direction_;

            //- Current power
            scalar power_;

            //- Index of the cell the ray currently occupies (or most recently occupied)
            //  It is -1 if it is not in a local processor cell
            label currentCell_;

            //- Number of bounces
            label bounceCount_;

            //- Global index of the ray
            label globalRayIndex_;

            //- Is it active?
            bool active_;

            //- List of all locations the ray has visited
            DynamicList<point> path_;


    // Public functions

    //- Runtime type information
    TypeName("compactRay");

    // Constructors

        //- Empty constructor
        compactRay();

        //- Construct from position, direction and power
        compactRay(const point& position, const vector& dir, scalar power);

        // Destructor

            virtual ~compactRay()
            {}


        // Public functions

        //- Reset the position, direction and power
        void reset(const point& position, const vector& dir, scalar power)
        {
            position_ = position;
            direction_ = dir / (mag(dir) + VSMALL);
            power_ = power;
            currentCell_ = -1;
            bounceCount_ = 0;
            globalRayIndex_ = -1;
            active_ = true;
        }

        //- Compare equals
        bool operator==(const compactRay& other) const
        {
            return
                position_ == other.position_
             && direction_ == other.direction_
             && power_ == other.power_
             && currentCell_ == other.currentCell_
             && bounceCount_ == other.bounceCount_
             && globalRayIndex_ == other.globalRayIndex_
             && active_ == other.active_
             && path_ == other.path_;
        }

        //- Compare not equals
        bool operator!=(const compactRay& other) const
        {
            return !(*this == other);
        }

        //- Output stream operator
        friend Ostream& operator<<(Ostream& os, const compactRay& ray)
        {
            return os  << ray.position_
                       << ray.direction_
                       << ray.power_
                       << ray.currentCell_
                       << ray.bounceCount_
                       << ray.globalRayIndex_
                       << ray.active_
                       << ray.path_;
        }

        //- Input stream operator
        friend Istream& operator>>(Istream& is, compactRay& ray)
        {
            return is  >> ray.position_
                       >> ray.direction_
                       >> ray.power_
                       >> ray.currentCell_
                       >> ray.bounceCount_
                       >> ray.globalRayIndex_
                       >> ray.active_
                       >> ray.path_;
        }
};


/*---------------------------------------------------------------------------*\
                        Helper functions and structs
\*---------------------------------------------------------------------------*/

struct combineRayLists
{
    void operator()
    (
        DynamicList<compactRay>& x,
        const DynamicList<compactRay>& y
    ) const
    {
        // Create a hash set of existing items for O(1) lookup
        HashSet<compactRay> existing;
        forAll(x, i)
        {
            existing.insert(x[i]);
        }

        // Add only new items
        forAll(y, i)
        {
            if (!existing.found(y[i]))
            {
                x.append(y[i]);
            }
        }
    }
};


struct combineLists
{
    void operator()
    (
        DynamicList<point>& x,
        const DynamicList<point>& y
    ) const
    {
        x.append(y);
    }
};


    struct combineRayPaths
{
    void operator()
    (
        DynamicList<DynamicList<point>>& x,
        const DynamicList<DynamicList<point>>& y
    ) const
    {
        forAll(y, i)
        {
            x.append(y[i]);
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
