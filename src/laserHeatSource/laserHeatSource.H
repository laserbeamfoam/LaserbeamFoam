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
    laserHeatSource

Description
    Calculate the deposition heat source field from a laser using ray-tracing.

    Laser parameters are defined via the provided dictionary.

    It is assumed that the laser enters the domain at the boundary with an
    outward-facing unit normal of (0 1 0).

    Laser rays are written as VTKs for visualisation.

Author
    Thomas Flint, The University of Manchester.
    Philip Cardiff, University College Dublin.

SourceFiles
    laserHeatSource.C

\*---------------------------------------------------------------------------*/

#ifndef laserHeatSource_H
#define laserHeatSource_H

#include "volFields.H"
#include "dictionary.H"
#include "interpolationTable.H"

#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class laserHeatSource Declaration
\*---------------------------------------------------------------------------*/


class CompactRay
{
    private:

    public:

    point origin_;
    vector direction_;
    scalar power_;
    scalar remainingPower_;
    label currentCell_;
    label bounceCount_;
    bool active_;
    
    CompactRay()
    : power_(0), remainingPower_(0), currentCell_(-1), 
      bounceCount_(0), active_(false) {}
      
    CompactRay(const point& origin, const vector& dir, scalar power)
    : origin_(origin), direction_(dir), power_(power), 
      remainingPower_(power), currentCell_(-1), bounceCount_(0), active_(true)
    {
        direction_ /= mag(direction_) + VSMALL;
    }
    void reset(const point& origin, const vector& dir, scalar power)
    {
        origin_ = origin;
        direction_ = dir / (mag(dir) + VSMALL);
        power_ = power;
        remainingPower_ = power;
        currentCell_ = -1;
        bounceCount_ = 0;
        active_ = true;
    }
};

class laserHeatSource
:
    public IOdictionary
{
    // Private data

        // Heat source field
        volScalarField deposition_;

        //- Field indicating the boundary faces where the laser enters the
        //  domain
        volScalarField laserBoundary_;

        //- Field for tracking error codes
        volScalarField errorTrack_;

        //- Field for tracking the latest ray to vist each cell
        volScalarField rayNumber_;

        //- Field for tracking the energy deposited in each cell
        volScalarField rayQ_;

        //- y dimension field
        volScalarField yDim_;

        //- Refine flag field
        volScalarField refineFlag_;

        //- Powder sim flag
        //  Needed for cases with spherical particles or spurious currents will
        //  smear out the particle interface un-physically
        const Switch powderSim_;

        //- new laser discretisation flag
        //  Tell the heat source if you want to discreetise the traditional laserbeamfoam way or in radial and polar coordinates
        const Switch Radial_Polar_HS_;

        //- List of the laser names
        wordList laserNames_;

        //- Copy of the laser dicts
        PtrList<dictionary> laserDicts_;

        //- List of the time vs laser position interpolation series
        PtrList<interpolationTable<vector>> timeVsLaserPosition_;

        //- List of the time vs laser power interpolation series
        PtrList<interpolationTable<scalar>> timeVsLaserPower_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        laserHeatSource(const laserHeatSource&);

        //- Disallow default bitwise assignment
        void operator=(const laserHeatSource&);

public:

    //- Runtime type information
    TypeName("laserHeatSource");

    // Constructors

        //- Construct from the mesh
        laserHeatSource
        (
            const fvMesh& mesh
        );

    // Destructor

        virtual ~laserHeatSource()
        {}


    // Member Functions

        //- Const access to the deposition field
        const volScalarField& deposition() const
        {
            return deposition_;
        }

        //- Non-const access to the yDim field
        volScalarField& yDim()
        {
            return yDim_;
        }

        //- Non-const access to the refineFlag field
        volScalarField& refineFlag()
        {
            return refineFlag_;
        }

        //- Const access to the powderSim flag
        const Switch& powderSim() const
        {
            return powderSim_;
        }

        //- Const access to the HS discretisation flag
        const Switch& Radial_Polar_HS() const
        {
            return Radial_Polar_HS_;
        }

        //- Update the deposition field for all lasers
        //  The filtered volume frame field (alphaFiltered) and the filtered
        //  interface normal field are required to define the interface
        //  The electrical resistivity field is passed as an argument, allowing
        //  a general laser heat source for multi-component and temperature
        //  dependent resistivity solvers
        //  This resets the deposition field to zero before updating it
        void updateDeposition
        (
            const volScalarField& alphaFiltered,
            const volVectorField& nFiltered,
            const volScalarField& resistivity_in
        );

        //- Update the deposition for the given laser parameters
        //  This adds to the existing deposition field and does not reset it to
        //  zero first
        void updateDeposition
        (
            const volScalarField& alphaFiltered,
            const volVectorField& nFiltered,
            const volScalarField& resistivity_in,
            const word& laserName,
            const vector& currentLaserPosition,
            const scalar currentLaserPower,
            const scalar laserRadius,
            const label N_sub_divisions,
            const label nRadial_,   //added for radial polar heat source
            const label nAngular_,   //added for radial polar heat source
            const vector& V_incident,
            const scalar wavelength,
            const scalar e_num_density,
            const scalar dep_cutoff,
            const scalar Radius_Flavour,
            const Switch useLocalSearch,
            const label maxLocalSearch
        );

        bool pointOnFace(const face& f, const point& p, const pointField& meshPoints);

        // void initialise_point_RP_linear(DynamicList<vector> *initial_points, DynamicList<scalar> *point_assoc_power);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
