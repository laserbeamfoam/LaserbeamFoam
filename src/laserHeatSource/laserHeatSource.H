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
#include "OFstream.H"
#include "compactRay.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                        Class laserHeatSource Declaration
\*---------------------------------------------------------------------------*/


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

        //- Powder sim flag
        //  Needed for cases with spherical particles or spurious currents will
        //  smear out the particle interface un-physically
        const Switch powderSim_;

        //- New laser discretisation flag
        //  Tell the heat source if you want to discretise the traditional
        //  laserbeamfoam way or in radial and polar coordinates
        const Switch radialPolarHeatSource_;

        //- List of the laser names
        wordList laserNames_;

        //- Copy of the laser dicts
        PtrList<dictionary> laserDicts_;

        //- List of the time vs laser position interpolation series
        PtrList<interpolationTable<vector>> timeVsLaserPosition_;

        //- List of the time vs laser power interpolation series
        PtrList<interpolationTable<scalar>> timeVsLaserPower_;

        //- Record the ray paths for each laser
        List<List<DynamicList<point>>> rayPaths_;

        //- Record of times when VTK files were written
        HashSet<scalar> vtkTimes_;

        //- CPU or GPU ray tracing
        enum rayTracingBackend
        {
            RT_CPU,
            RT_GPU,
            RT_AUTO
        };
        rayTracingBackend rayBackend_;

        //- Is a GPU available?
        bool gpuAvailable_;

        //- Ray tracing algorithm
        enum cellSearchMode
        {
            CSM_LEGACY,
            CSM_VOXEL,
            CSM_VOXEL_WITH_FALLBACK
        };
        cellSearchMode cellSearchMode_;

        //- Has voxel data been built yet?
        bool voxelInitialised_;

        //- Voxel-based ray tracing parameter: number of voxels in each direction
        label nVoxelsX_, nVoxelsY_, nVoxelsZ_;

        //- Voxel-based ray tracing parameter: Voxel edge coordinates
        scalarField xEdges_, yEdges_, zEdges_;

        //- Voxel-based ray tracing parameter: Flattened voxel
        //  Cell map: size = nVoxelsX_*nVoxelsY_*nVoxelsZ_
        labelList voxelCell_;

        //- Voxel-based ray tracing parameter
        //  Per-cell axis-aligned bounding boxes
        List<vector> cellMin_;
        List<vector> cellMax_;

        //- Global bounding box of the mesh
        boundBox globalBB_;


    // Private Member Functions

        //- Create the initial rays
        //  This involves creating a list of points (starting position of the
        //  laser rays) and their associated power (J/s)
        //  The rays list should be the same on all
        //  processors when returned from this function
        void createInitialRays
        (
            List<compactRay>& rays,
            const fvMesh& mesh,
            const vector& currentLaserPosition,
            const scalar laserRadius,
            const label N_sub_divisions,
            const label nRadial,
            const label nAngular,
            const vector& V_incident,
            const scalar Radius_Flavour,
            const scalar Q_cond,
            const scalar beam_radius
        ) const;

        //- Create voxel data for voxel-based ray tracing
        void initialiseVoxelData
        (
            const boundBox& globalBB,
            const label Nx,
            const label Ny,
            const label Nz,
            const label maxLocalSearch
        );

        //- Wrapper for ray cell search
        label findCellForRay
        (
            const point& p,
            const label currentCell,
            const fvMesh& mesh,
            const label maxLocalSearch,
            const bool debug = false
        ) const;

        //- Voxel-based search
        label findCellVoxelised
        (
            const point& p,
            const label currentCell,
            const fvMesh& mesh,
            const label maxLocalSearch,
            const bool debug = false
        ) const;

        //- Helper: map (i,j,k) to flat index
        inline label voxelIndex
        (
            const label i, const label j, const label k
        ) const
        {
            return (k*nVoxelsY_ + j)*nVoxelsX_ + i;
        }

        //- Helper: compute voxel indices from a point
        inline bool voxelIndices
        (
            const point& p, label& i, label& j, label& k
        ) const;

        //- Detech if a GPU is available
        bool detectGPU();

        //- Propagate rays on CPU
        void propagateRaysCPU
        (
            DynamicList<compactRay>& remainingGlobalRays,
            const fvMesh& mesh,
            const scalarField& VI,
            const volScalarField& alphaFiltered,
            const volVectorField& nFiltered,
            const volScalarField& resistivity_in,
            const scalar plasma_frequency,
            const scalar angular_frequency,
            const scalar dep_cutoff,
            const scalar rayPowerAbsTol,
            const boundBox& globalBB,
            const label maxLocalSearch,
            const label laserID
        );

        //- Propagate rays on CPU
        void propagateRaysGPU
        (
            DynamicList<compactRay>& remainingGlobalRays,
            const fvMesh& mesh,
            const scalarField& VI,
            const volScalarField& alphaFiltered,
            const volVectorField& nFiltered,
            const volScalarField& resistivity_in,
            const scalar plasma_frequency,
            const scalar angular_frequency,
            const scalar dep_cutoff,
            const scalar rayPowerAbsTol,
            const boundBox& globalBB,
            const label maxLocalSearch,
            const label laserID
        );

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

        //- Const access to the powderSim flag
        const Switch& powderSim() const
        {
            return powderSim_;
        }

        //- Const access to the HS discretisation flag
        const Switch& radialPolarHeatSource() const
        {
            return radialPolarHeatSource_;
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
            const label laserID,
            const vector& currentLaserPosition,
            const scalar currentLaserPower,
            const scalar laserRadius,
            const label N_sub_divisions,
            const label nRadial,
            const label nAngular,
            const vector& V_incident,
            const scalar wavelength,
            const scalar e_num_density,
            const scalar dep_cutoff,
            const scalar Radius_Flavour,
            const Switch useLocalSearch,
            const label maxLocalSearch,
            const scalar rayPowerRelTol,
            const boundBox& globalBB
        );

        //- Update the deposition for the given laser parameters using the CPU
        //  implementation
        void updateDepositionCPU
        (
            const volScalarField& alphaFiltered,
            const volVectorField& nFiltered,
            const volScalarField& resistivity_in,
            const label laserID,
            const vector& currentLaserPosition,
            const scalar currentLaserPower,
            const scalar laserRadius,
            const label N_sub_divisions,
            const label nRadial,
            const label nAngular,
            const vector& V_incident,
            const scalar wavelength,
            const scalar e_num_density,
            const scalar dep_cutoff,
            const scalar Radius_Flavour,
            const Switch useLocalSearch,
            const label maxLocalSearch,
            const scalar rayPowerRelTol,
            const boundBox& globalBB
        );

        //- Update the deposition for the given laser parameters using the GPU
        //  implementation
        void updateDepositionGPU
        (
            const volScalarField& alphaFiltered,
            const volVectorField& nFiltered,
            const volScalarField& resistivity_in,
            const label laserID,
            const vector& currentLaserPosition,
            const scalar currentLaserPower,
            const scalar laserRadius,
            const label N_sub_divisions,
            const label nRadial,
            const label nAngular,
            const vector& V_incident,
            const scalar wavelength,
            const scalar e_num_density,
            const scalar dep_cutoff,
            const scalar Radius_Flavour,
            const Switch useLocalSearch,
            const label maxLocalSearch,
            const scalar rayPowerRelTol,
            const boundBox& globalBB
        );

        //- Const access to the global bounding box
        const boundBox& globalBB() const
        {
            return globalBB_;
        }

        //- Write the given ray paths to a VTK file with the given name
        void writeRayPathsToVTK
        (
            const List<DynamicList<point>>& rays,
            const fileName& filename
        );

        //- Write rayPath_ to a VTK file with the given name
        void writeRayPathsToVTK();

        //- Write a "vtk.series" file with details about the VTK files written
        //  so far. This file can be conveniently opened in ParaView and will
        //  load the rays at the correct physical times
        void writeRayPathVTKSeriesFile() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
