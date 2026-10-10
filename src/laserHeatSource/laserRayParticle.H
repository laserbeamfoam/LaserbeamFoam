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
    Foam::laserRayParticle

Description
    Lagrangian particle representing a single laser ray for use in a
    Discrete Transfer Radiation Method (DTRM) style ray-tracing algorithm.

    Each particle carries a ray direction, power, and area element.
    The particle uses OpenFOAM's built-in particle tracking to traverse
    cells exactly (face-to-face), with automatic parallel transfer at
    processor boundaries handled by the Cloud infrastructure.

    Inspired by OpenFOAM's DTRMParticle used in the laserDTRM radiation
    model.

Author
    Tom Flint, The University of Manchester.
    Philip Cardiff, University College Dublin.

SourceFiles
    laserRayParticle.C

\*---------------------------------------------------------------------------*/

#ifndef laserRayParticle_H
#define laserRayParticle_H

#include "particle.H"
#include "Cloud.H"
#include "volFields.H"
#include "DynamicList.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
class laserRayParticle;

// Cloud typedef
typedef Cloud<laserRayParticle> laserRayParticleCloud;

// Friend operators
Ostream& operator<<(Ostream&, const laserRayParticle&);


/*---------------------------------------------------------------------------*\
                      Class laserRayParticle Declaration
\*---------------------------------------------------------------------------*/

class laserRayParticle
:
    public particle
{

public:

    // Public data carried by the ray

        //- Current ray direction (unit vector)
        vector direction_;

        //- Current ray power [W]
        scalar power_;

        //- Initial ray power [W] (for tolerance checks)
        scalar initialPower_;

        //- Area element associated with this ray [m^2]
        scalar dA_;

        //- Number of reflections
        label bounceCount_;

        //- Smallest absorptivity encountered by this ray
        scalar minAbsorptivity_;

        //- Number of near-zero absorptivity reflections seen by this ray
        label nearZeroAbsorptivityCount_;

        //- Global ray index (for VTK output ordering)
        label globalRayIndex_;

        //- Is the ray still active?
        bool active_;

        //- Accumulated path points for VTK visualisation
        DynamicList<point> path_;


    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
    public:

        // References to fields needed during tracking

            //- Filtered volume fraction field
            const volScalarField& alphaFiltered_;

            //- Filtered interface normal field
            const volVectorField& nFiltered_;

            //- Electrical resistivity field
            const volScalarField& resistivity_;

            //- Reference to the deposition field (accumulated into)
            volScalarField& deposition_;

            //- Reference to the ray power visualisation field
            volScalarField& rayQ_;

            //- Reference to the ray number visualisation field
            volScalarField& rayNumber_;

            //- Deposition alpha cutoff
            scalar depCutoff_;

            //- Plasma frequency [rad/s]
            scalar plasmaFrequency_;

            //- Angular frequency of laser light [rad/s]
            scalar angularFrequency_;

            //- Ray power tolerance relative to each ray's initial power
            scalar rayPowerRelTol_;

            //- Maximum number of reflections before handling a ray
            label maxRayBounces_;

            //- Action when maxRayBounces is exceeded
            word maxRayBounceAction_;

            //- Ray absorptivity below this value is counted as near-zero
            scalar nearZeroAbsorptivityTol_;

            //- Number of rays handled by maxRayBounces on this rank
            label& maxBounceTerminations_;

            //- Residual ray power handled by maxRayBounces on this rank
            scalar& maxBounceResidualPower_;

            //- Smallest absorptivity encountered on this rank
            scalar& minAbsorptivity_;

            //- Number of near-zero absorptivity interactions on this rank
            label& nearZeroAbsorptivityCount_;

            //- Largest bounce count encountered on this rank
            label& maxObservedBounces_;

            //- Collected ray paths from finished particles
            //  (particles store their paths here before dying,
            //   since dead particles are removed from the cloud)
            DynamicList<label>& finishedRayIDs_;
            DynamicList<DynamicList<point>>& finishedRayPaths_;


        // Constructor
        trackingData
        (
            const Cloud<laserRayParticle>& cloud,
            const volScalarField& alphaFiltered,
            const volVectorField& nFiltered,
            const volScalarField& resistivity,
            volScalarField& deposition,
            volScalarField& rayQ,
            volScalarField& rayNumber,
            const scalar depCutoff,
            const scalar plasmaFrequency,
            const scalar angularFrequency,
            const scalar rayPowerRelTol,
            const label maxRayBounces,
            const word& maxRayBounceAction,
            const scalar nearZeroAbsorptivityTol,
            label& maxBounceTerminations,
            scalar& maxBounceResidualPower,
            scalar& minAbsorptivity,
            label& nearZeroAbsorptivityCount,
            label& maxObservedBounces,
            DynamicList<label>& finishedRayIDs,
            DynamicList<DynamicList<point>>& finishedRayPaths
        )
        :
            particle::trackingData(cloud),
            alphaFiltered_(alphaFiltered),
            nFiltered_(nFiltered),
            resistivity_(resistivity),
            deposition_(deposition),
            rayQ_(rayQ),
            rayNumber_(rayNumber),
            depCutoff_(depCutoff),
            plasmaFrequency_(plasmaFrequency),
            angularFrequency_(angularFrequency),
            rayPowerRelTol_(rayPowerRelTol),
            maxRayBounces_(maxRayBounces),
            maxRayBounceAction_(maxRayBounceAction),
            nearZeroAbsorptivityTol_(nearZeroAbsorptivityTol),
            maxBounceTerminations_(maxBounceTerminations),
            maxBounceResidualPower_(maxBounceResidualPower),
            minAbsorptivity_(minAbsorptivity),
            nearZeroAbsorptivityCount_(nearZeroAbsorptivityCount),
            maxObservedBounces_(maxObservedBounces),
            finishedRayIDs_(finishedRayIDs),
            finishedRayPaths_(finishedRayPaths)
        {}
    };


    // Static data members

        //- Runtime type information
        TypeName("laserRayParticle");


    // Constructors

        //- Construct from components
        laserRayParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const label cellI,
            const vector& direction,
            const scalar power,
            const scalar dA,
            const label globalRayIndex
        );

        //- Construct from Istream
        laserRayParticle
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct and return a clone
        autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new laserRayParticle(*this));
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<laserRayParticle> operator()(Istream& is) const
            {
                return autoPtr<laserRayParticle>
                (
                    new laserRayParticle(mesh_, is, true)
                );
            }
        };


    //- Destructor
    virtual ~laserRayParticle() = default;


    // Member Functions

        //- Move the particle (called by Cloud::move)
        //  This is the main ray-tracing logic: track face-to-face,
        //  check for interfaces, deposit energy, reflect.
        bool move(Cloud<laserRayParticle>& cloud, trackingData& td, const scalar trackTime);

        //- Overridable function to handle the particle hitting a patch
        //  Return true if handled (skip specific patch handler), false
        //  to continue to hitProcessorPatch/hitWallPatch etc.
        bool hitPatch(Cloud<laserRayParticle>& cloud, trackingData& td);

        //- Hit a processor patch - transfer to the neighbour processor
        void hitProcessorPatch
        (
            Cloud<laserRayParticle>& cloud,
            trackingData& td
        );

        //- Hit a wall patch - absorb remaining power
        void hitWallPatch
        (
            Cloud<laserRayParticle>& cloud,
            trackingData& td
        );

        //- Compute Fresnel absorptivity at an interface
        static scalar computeFresnelAbsorptivity
        (
            const vector& rayDirection,
            const vector& surfaceNormal,
            const scalar plasmaFrequency,
            const scalar angularFrequency,
            const scalar resistivity
        );

        //- Read/write fields for parallel transfer
        static void readFields(Cloud<laserRayParticle>& cloud);
        static void writeFields(const Cloud<laserRayParticle>& cloud);


    // Ostream operator

        friend Ostream& operator<<
        (
            Ostream& os,
            const laserRayParticle& p
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
