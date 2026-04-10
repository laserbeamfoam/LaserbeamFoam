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

\*---------------------------------------------------------------------------*/

#include "laserHeatSource.H"
#include "fvc.H"
#include "constants.H"
#include "SortableList.H"

#include "findLocalCell.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(laserHeatSource, 0);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


void laserHeatSource::seedRayCloud
(
    Cloud<laserRayParticle>& cloud,
    const fvMesh& mesh,
    const vector& currentLaserPosition,
    const scalar laserRadius,
    const label N_sub_divisions,
    const label nRadial,
    const label nAngular,
    const vector& V_incident,
    const scalar Radius_Flavour,
    const scalar Q_cond,
    const scalar beam_radius,
    label& nTotalRays
) const
{
    // ================================================================
    // Compute ALL ray starting positions and powers identically on
    // every processor (same as the old createInitialRays), but only
    // add a particle to the cloud on the processor that owns the
    // starting cell.
    // ================================================================

    const scalarField& yDimI = yDim_;
    const scalar pi = constant::mathematical::pi;

    // -- Build global list of ray positions and powers --

    DynamicList<vector> initial_points;
    DynamicList<scalar> point_assoc_power;

    // All seed points will be perturbed by a small amount to avoid landing
    // exactly on a face
    // The magnitude of the perturbation is calculated as a small fraction of
    // the smallest cell dimension, as opposed to an absolute value, and we take
    // care to perturb in the ray direction and transverse to the ray direction

    // Ray direction
    const vector d = V_incident/(mag(V_incident) + VSMALL);

    // Pick a stable “arbitrary” axis not parallel to d
    const vector a =
        (mag(d & vector(1,0,0)) < 0.9) ? vector(1,0,0) : vector(0,0,1);

    // Perpendicular unit vector
    vector t = d ^ a;
    t /= (mag(t) + VSMALL);

    // Minimum length scale
    const scalar lCell = cbrt(gMin(mesh.V()));
    const scalar eps = sqrt(SMALL)*lCell;

    // Safe perturbation in the ray direction and transverse ray direction
    const vector perturbation = eps*(d + t);
    Info<< "    Perturbation vector for ray seed points = " << perturbation
        << endl;

    if (radialPolarHeatSource())
    {
        Info<< "nRadial: " << nRadial << nl
            << "nAngular: " << nAngular << endl;

        const scalar rMax = 1.5*beam_radius;

        // ================================================================
        // Improved radial discretisation for smoother Gaussian
        //
        // Use cell-edge based radial boundaries to compute exact annulus
        // areas. Each ring iR covers r from r_inner[iR] to r_outer[iR].
        // The sample point is at the area-weighted centroid of the annulus.
        // ================================================================

        // Compute radial boundaries (edges of annular rings)
        List<scalar> radialBoundaries(nRadial + 1);
        for (label iR = 0; iR <= nRadial; ++iR)
        {
            // Uniform spacing of boundaries
            radialBoundaries[iR] = rMax * scalar(iR) / scalar(nRadial);
        }

        // Compute sample points (area-weighted centroid of each annulus)
        // For an annulus from r1 to r2, the centroid is at:
        //   r_centroid = (2/3) * (r2^3 - r1^3) / (r2^2 - r1^2)
        List<scalar> radialPoints(nRadial);
        List<scalar> annulusAreas(nRadial);

        for (label iR = 0; iR < nRadial; ++iR)
        {
            const scalar r_inner = radialBoundaries[iR];
            const scalar r_outer = radialBoundaries[iR + 1];

            // Exact annulus area (full ring, will be divided by nAngular)
            annulusAreas[iR] = pi * (sqr(r_outer) - sqr(r_inner));

            // Area-weighted centroid for sample position
            if (r_inner < SMALL)
            {
                // Central disc: centroid at 2/3 * r_outer
                radialPoints[iR] = (2.0/3.0) * r_outer;
            }
            else
            {
                // Annulus: use exact centroid formula
                radialPoints[iR] =
                    (2.0/3.0)
                   * (pow3(r_outer) - pow3(r_inner))
                   / (sqr(r_outer) - sqr(r_inner));
            }
        }

        const label totalSamples = nRadial * nAngular;
        const label samplesPerProc = totalSamples/Pstream::nProcs();
        const label remainder = totalSamples % Pstream::nProcs();
        const label myRank = Pstream::myProcNo();
        const label startIdx =
            myRank * samplesPerProc + min(myRank, remainder);
        const label endIdx =
            startIdx + samplesPerProc + (myRank < remainder ? 1 : 0);
        const label localSamples = endIdx - startIdx;

        const point P0
        (
            currentLaserPosition.x(),
            currentLaserPosition.y(),
            currentLaserPosition.z()
        );

        const vector V_i(V_incident/(mag(V_incident) + SMALL));

        const vector a =
            (mag(V_i.z()) < 0.9) ? vector(0, 0, 1) : vector(0, 1, 0);
        vector u = (V_i ^ a);
        u = u/mag(u);
        const vector v = (V_i ^ u);

        for (label localIdx = 0; localIdx < localSamples; ++localIdx)
        {
            const label globalIdx = startIdx + localIdx;
            const label iTheta = globalIdx/nRadial;
            const label iR = globalIdx % nRadial;

            const scalar theta = 2.0*pi*iTheta/nAngular;
            const scalar r = radialPoints[iR];

            // Area element for this ray (annulus sector)
            const scalar area = annulusAreas[iR] / scalar(nAngular);

            const scalar x_local = r*cos(theta);
            const scalar y_local = r*sin(theta);

            const vector globalPos = P0 + x_local*u + y_local*v;

            initial_points.append(globalPos + perturbation);

            // Gaussian power distribution using the sample point radius
            point_assoc_power.append
            (
                area
               *(
                    Radius_Flavour*Q_cond
                   /(Foam::pow(beam_radius, 2.0)*pi)
                )
               *Foam::exp
                (
                  - Radius_Flavour
                   *(Foam::pow(r, 2.0)/Foam::pow(beam_radius, 2.0))
                )
            );
        }
    }
    else // One ray for each boundary patch face within the laser radius
    {
        const vectorField& CI = mesh.C();

        forAll(CI, celli)
        {
            const scalar x_coord = CI[celli].x();
            const scalar z_coord = CI[celli].z();

            const scalar r =
                sqrt
                (
                    sqr(x_coord - currentLaserPosition.x())
                  + sqr(z_coord - currentLaserPosition.z())
                );

            if
            (
                r <= (1.5*beam_radius)
             && laserBoundary_[celli] > SMALL
            )
            {
                for (label Ray_j = 0; Ray_j < N_sub_divisions; Ray_j++)
                {
                    for (label Ray_k = 0; Ray_k < N_sub_divisions; Ray_k++)
                    {
                        point p_1
                        (
                            CI[celli].x()
                          - (yDimI[celli]/2.0)
                          + ((yDimI[celli]/(N_sub_divisions+1))*(Ray_j+1)),
                            CI[celli].y(),
                            CI[celli].z()
                          - (yDimI[celli]/2.0)
                          + ((yDimI[celli]/(N_sub_divisions+1))*(Ray_k+1))
                        );

                        // Perturb the point
                        p_1 += perturbation;

                        initial_points.append(p_1);

                        point_assoc_power.append
                        (
                            sqr(yDimI[celli]/N_sub_divisions)
                           *(
                                (Radius_Flavour*Q_cond)
                               /(Foam::pow(beam_radius, 2.0)*pi)
                            )
                           *Foam::exp
                            (
                              - Radius_Flavour
                               *(
                                    Foam::pow(r, 2.0)
                                   /Foam::pow(beam_radius, 2.0)
                                )
                            )
                        );
                    }
                }
            }
        }
    }

    // List with size equal to number of processors
    // Gather all initial points/powers to all processors
    List<pointField> gatheredData(Pstream::nProcs());
    List<scalarField> gatheredPowers(Pstream::nProcs());

    // Populate and gather the list onto the master processor.
    // Distibulte the data accross the different processors
    gatheredData[Pstream::myProcNo()] = initial_points;
    Pstream::gatherList(gatheredData);
    Pstream::broadcastList(gatheredData);

    gatheredPowers[Pstream::myProcNo()] = point_assoc_power;
    Pstream::gatherList(gatheredPowers);
    Pstream::broadcastList(gatheredPowers);

    // List of initial points
    pointField rayCoords
    (
        ListListOps::combine<Field<vector>>
        (
            gatheredData,
            accessOp<Field<vector>>()
        )
    );

    scalarField rayPowers
    (
        ListListOps::combine<Field<scalar>>
        (
            gatheredPowers,
            accessOp<Field<scalar>>()
        )
    );


    // Create a list of Ray objects
    nTotalRays = rayCoords.size();

    label seedCellI = -1;

    forAll(rayCoords, i)
    {
        const label cellI = findLocalCell
            (
                rayCoords[i],
                seedCellI,
                mesh,
                100,    // maxLocalSearch
                false   // debug
            );

        if (cellI >= 0)
        {
            // Update seed for next iteration - nearby rays benefit
            seedCellI = cellI;

            // This processor owns this ray's starting cell
            laserRayParticle* pPtr = new laserRayParticle
                (
                    mesh,
                    rayCoords[i],
                    cellI,
                    V_incident,
                    rayPowers[i],
                    0.0,        // dA
                    i           // globalRayIndex
                );

            cloud.addParticle(pPtr);
        }
    }

    Info<< "    Total rays: " << nTotalRays
        << ", local particles: " << cloud.size() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

laserHeatSource::laserHeatSource
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "LaserProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    deposition_
    (
        IOobject
        (
            "Deposition",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("deposition", dimensionSet(1, -1, -3, -0, 0), -1.0)
    ),
    laserBoundary_
    (
        IOobject
        (
            "Laser_boundary",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    errorTrack_
    (
        IOobject
        (
            "errorTrack",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "errorTrack", dimensionSet(0, 0, 0, -0, 0), 0.0
        )
    ),
    rayNumber_
    (
        IOobject
        (
            "rayNumber",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rayNumber", dimensionSet(0, 0, 0, -0, 0), -1.0)
    ),
    rayQ_
    (
        IOobject
        (
            "rayQ",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("rayQ", dimensionSet(1, 0, -3, 0, 0), scalar(0.0))
    ),
    yDim_
    (
        IOobject
        (
            "yDim",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("yDim", dimensionSet(0, 1, 0, 0, 0), 1.0)
    ),
    powderSim_(lookupOrDefault<Switch>("PowderSim", false)), // Needed for cases where we have spherical particles (that have jagged edges or cartesian grids) which will cause spurious currents due to these large gradients - we extend the damping term by 1 cell
    radialPolarHeatSource_
    (
        found("radialPolarHeatSource")
      ? Switch(lookup("radialPolarHeatSource"))
      : lookupOrDefault<Switch>("Radial_Polar_HS", true)
    ),
    laserNames_(0),
    laserDicts_(0),
    timeVsLaserPosition_(0),
    timeVsLaserPower_(0),
    rayPaths_(0),
    vtkTimes_(),
    globalBB_(mesh.bounds())
{
    Info<< "radialPolarHeatSource = " << radialPolarHeatSource_ << endl;

    // Calculate global bounding box
    {
        boundBox localBB = mesh.bounds();
        globalBB_ = localBB;

        reduce(globalBB_.min(), minOp<vector>());
        reduce(globalBB_.max(), maxOp<vector>());

        globalBB_.inflate(0.01);
        // Inflate the bounding box by 1% to avoid issues with rays starting just at the mesh boundary -> probably not an issue but worth doing

        Info<< "Scaled global mesh bounding box: " << globalBB_ << endl;
    }

    // Initialise the laser power and position
    if (found("lasers"))
    {
        const PtrList<entry> laserEntries(lookup("lasers"));

        laserNames_.setSize(laserEntries.size());
        laserDicts_.setSize(laserEntries.size());
        timeVsLaserPosition_.setSize(laserEntries.size());
        timeVsLaserPower_.setSize(laserEntries.size());

        if (Pstream::master())
        {
            rayPaths_.setSize(laserEntries.size());
        }

        forAll(laserEntries, laserI)
        {
            laserNames_[laserI] = laserEntries[laserI].keyword();
            Info<< "Reading laser " << laserNames_[laserI] << endl;

            laserDicts_.set
            (
                laserI, new dictionary(laserEntries[laserI].dict())
            );

            timeVsLaserPosition_.set
            (
                laserI,
                new interpolationTable<vector>
                (
                    laserEntries[laserI].dict().subDict("timeVsLaserPosition")
                )
            );

            timeVsLaserPower_.set
            (
                laserI,
                new interpolationTable<scalar>
                (
                    laserEntries[laserI].dict().subDict("timeVsLaserPower")
                )
            );
        }

            // Check that a single laser is not also defined

        if (found("timeVsLaserPosition"))
        {
            FatalErrorInFunction
                << "timeVsLaserPosition should not be defined in the main dict"
                << " if a list of lasers is provided" << exit(FatalError);
        }

        if (found("timeVsLaserPower"))
        {
            FatalErrorInFunction
                << "timeVsLaserPower should not be defined in the main dict"
                << " if a list of lasers is provided" << exit(FatalError);
        }
    }
    else
    {
        // There is no lists of lasers, just one
        // Backward compatibility: single laser specified in main dict
        rayPaths_.setSize(1);
        laserNames_.setSize(1);
        laserDicts_.setSize(1);
        timeVsLaserPosition_.setSize(1);
        timeVsLaserPower_.setSize(1);

        laserNames_[0] = "laser0";

        // Copy the main dict
        laserDicts_.set(0, new dictionary(*this));

        timeVsLaserPosition_.set
        (
            0,
            new interpolationTable<vector>(subDict("timeVsLaserPosition"))
        );

        timeVsLaserPower_.set
        (
            0,
            new interpolationTable<scalar>(subDict("timeVsLaserPower"))
        );
    }

        // Update laserBoundary - only used in old boundary based initialisation
    laserBoundary_ = fvc::average(laserBoundary_);

    if (debug)
    {
        errorTrack_.writeOpt() = IOobject::AUTO_WRITE;
        rayNumber_.writeOpt() = IOobject::AUTO_WRITE;
    }

    // Give errors if the old input format is found
    if (found("HS_bg"))
    {
        FatalErrorInFunction
            << "'HS_bg' is deprecated: please instead specify the laser "
            << "position in time via the laserPositionVsTime sub-dict"
            << exit(FatalError);
    }
    // Give errors if the old input format is found
    if (found("HS_lg"))
    {
        FatalErrorInFunction
            << "'HS_lg' is deprecated: please instead specify the laser "
            << "position in time via the laserPositionVsTime sub-dict"
            << exit(FatalError);
    }
    // Give errors if the old input format is found
    if (found("HS_velocity"))
    {
        FatalErrorInFunction
            << "'HS_velocity' is deprecated: please instead specify the laser "
            << "position in time via the laserPositionVsTime sub-dict"
            << exit(FatalError);
    }
    // Give errors if the old input format is found
    if (found("HS_Q"))
    {
        FatalErrorInFunction
            << "'HS_Q' is deprecated: please instead specify the laser "
            << "power in time via the laserPowereVsTime sub-dict"
            << exit(FatalError);
    }
    // Give errors if the old input format is found
    if (found("elec_resistivity"))
    {
        FatalErrorInFunction
            << "'elec_resistivity' is deprecated: resistivity is now "
            << "passed in from the solver as a field"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void laserHeatSource::updateDeposition
(
    const volScalarField& alphaFiltered,
    const volVectorField& nFiltered,
    const volScalarField& resistivity_in
)
{
    Info<< "Updating deposition" << endl;

    deposition_ *= 0.0;
    laserBoundary_ *= 0.0;
    laserBoundary_ = fvc::average(laserBoundary_);
    errorTrack_ *= 0.0;
    rayNumber_ *= 0.0;
    rayQ_ *= 0.0;

    const scalar time = deposition_.time().value();

    forAll(laserNames_, laserI)
    {
        // Lookup the current laser position and power
        vector currentLaserPosition =
            timeVsLaserPosition_[laserI](time);
        const scalar currentLaserPower =
            timeVsLaserPower_[laserI](time);

        Info<< "Laser: " << laserNames_[laserI] << nl
            << "    mean position = " << currentLaserPosition << nl
            << "    power = " << currentLaserPower << endl;

            // Dict for current laser
        const dictionary& dict = laserDicts_[laserI];

        // If defined, add oscillation to laser position
        if (dict.found("HS_oscAmpX"))
        {
            const scalar oscAmpX(readScalar(dict.lookup("HS_oscAmpX")));
            const scalar oscFreqX(readScalar(dict.lookup("HS_oscFreqX")));
            const scalar pi = constant::mathematical::pi;
            currentLaserPosition[vector::X] +=
                oscAmpX*sin(2*pi*oscFreqX*time);
        }// If defined, add oscillation to laser position

        if (dict.found("HS_oscAmpZ"))
        {
            const scalar oscAmpZ(readScalar(dict.lookup("HS_oscAmpZ")));
            const scalar oscFreqZ(readScalar(dict.lookup("HS_oscFreqZ")));
            const scalar pi = constant::mathematical::pi;
            currentLaserPosition[vector::Z] +=
                oscAmpZ*cos(2*pi*oscFreqZ*time);
        }// If defined, add oscillation to laser position

        Info<< "    position including any oscillation = "
            << currentLaserPosition << endl;

        scalar laserRadius = 0.0;
        if (dict.found("HS_a") && dict.found("laserRadius"))
        {
            FatalErrorInFunction
                << "The laser radius should be specified via 'laserRadius' or "
                << "'HS_a', not both!" << exit(FatalError);
        }

        if (dict.found("HS_a"))
        {
            laserRadius = readScalar(dict.lookup("HS_a"));
        }
        else if (dict.found("laserRadius"))
        {
            laserRadius = readScalar(dict.lookup("laserRadius"));
        }
        else
        {
            FatalErrorInFunction
                << "The laser radius should be specified via 'laserRadius' "
                << "or 'HS_a'"
                << exit(FatalError);
        }

        const label nRadial(dict.lookupOrDefault<label>("nRadial", 5));
        const label nAngular(dict.lookupOrDefault<label>("nAngular", 30));
        const label N_sub_divisions
        (
            dict.lookupOrDefault<label>("N_sub_divisions", 1)
        );
        const vector V_incident(dict.lookup("V_incident"));
        const scalar wavelength(readScalar(dict.lookup("wavelength")));
        const scalar e_num_density
        (
            readScalar(dict.lookup("e_num_density"))
        );
        const scalar dep_cutoff
        (
            dict.lookupOrDefault<scalar>("dep_cutoff", 0.5)
        );
        const scalar Radius_Flavour
        (
            dict.lookupOrDefault<scalar>("Radius_Flavour", 2.0)
        );
        const Switch useLocalSearch
        (
            dict.lookupOrDefault<Switch>("useLocalSearch", true)
        );
        const label maxLocalSearch
        (
            dict.lookupOrDefault<label>("maxLocalSearch", 100)
        );
        const scalar rayPowerRelTol
        (
            dict.lookupOrDefault<scalar>("rayPowerRelTol", 1e-6)
        );

        updateDeposition
        (
            alphaFiltered,
            nFiltered,
            resistivity_in,
            laserI,
            currentLaserPosition,
            currentLaserPower,
            laserRadius,
            N_sub_divisions,
            nRadial,
            nAngular,
            V_incident,
            wavelength,
            e_num_density,
            dep_cutoff,
            Radius_Flavour,
            useLocalSearch,
            maxLocalSearch,
            rayPowerRelTol,
            globalBB_
        );
    }

    deposition_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Lagrangian particle-based parallel ray-tracing: updateDeposition
//
//  Algorithm:
//    1. Create a Cloud<laserRayParticle>
//    2. Seed it with ray particles (one per ray, owned by local proc)
//    3. Call cloud.move() - OpenFOAM handles:
//       - Exact face-to-face tracking (no step size, no missed cells)
//       - Automatic transfer at processor boundaries
//       - Hit callbacks for walls, patches, processor patches
//    4. Collect ray paths from finished particles for VTK output
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void laserHeatSource::updateDeposition
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
)
{
    const fvMesh& mesh = deposition_.mesh();
    const scalar pi = constant::mathematical::pi;
    const scalar beam_radius = laserRadius;

    // Material constants
    const scalar plasma_frequency = Foam::sqrt
    (
        (
            e_num_density
           *constant::electromagnetic::e.value()
           *constant::electromagnetic::e.value()
        )
       /(
           constant::atomic::me.value()
          *constant::electromagnetic::epsilon0.value()
       )
    );
    const scalar angular_frequency =
        2.0*pi*constant::universal::c.value()/wavelength;

    // ==================================================================
    // Step 1: Create the ray particle cloud
    // ==================================================================

    Cloud<laserRayParticle> rayCloud
    (
        mesh,
        "laserRays",
        IDLList<laserRayParticle>()
    );

    // ==================================================================
    // Step 2: Seed the cloud with ray particles
    //
    // seedRayCloud computes ALL ray positions/powers globally, but
    // each processor only adds particles for cells it owns.
    // ==================================================================

    label nTotalRays = 0;

    seedRayCloud
    (
        rayCloud,
        mesh,
        currentLaserPosition,
        laserRadius,
        N_sub_divisions,
        nRadial,
        nAngular,
        V_incident,
        Radius_Flavour,
        currentLaserPower,
        beam_radius,
        nTotalRays
    );

    // Compute absolute power tolerance from relative tolerance
    // (Based on average ray power, not total power)
    const scalar rayPowerAbsTol =
        rayPowerRelTol*currentLaserPower/max(nTotalRays, 1);

    // ==================================================================
    // Step 3: Move the particles through the mesh
    //
    //  Cloud::move() handles everything:
    //    - Face-to-face tracking via trackToAndHitFace()
    //    - Automatic parallel transfer at processor patches
    //    - Callback to hitWallPatch() / hitPatch() at boundaries
    //
    //  The deposition field is written to directly by each particle
    //  during tracking - only to local cells, so no reduce needed.
    // ==================================================================

    DynamicList<label> finishedRayIDs;
    DynamicList<DynamicList<point>> finishedRayPaths;

    laserRayParticle::trackingData td
    (
        rayCloud,
        alphaFiltered,
        nFiltered,
        resistivity_in,
        deposition_,
        rayQ_,
        rayNumber_,
        dep_cutoff,
        plasma_frequency,
        angular_frequency,
        rayPowerAbsTol,
        finishedRayIDs,
        finishedRayPaths
    );

    // Move all ray particles through the mesh.
    // trackTime is not physically meaningful for rays (instantaneous),
    // but the Cloud::move interface requires it. We use GREAT to ensure
    // the rays are tracked to completion in a single call.
    rayCloud.storeGlobalPositions();
    rayCloud.move(rayCloud, td, GREAT);

    // ==================================================================
    // Step 4: Gather ray path segments to master for VTK output
    //
    // Each ray may produce multiple segments (one per processor it
    // traverses). Segments are stored in trackingData both at
    // processor-boundary crossings and at particle death. We
    // gather all segments and concatenate them per rayID.
    // ==================================================================

    if (Pstream::master())
    {
        rayPaths_[laserID].clear();
        rayPaths_[laserID].setSize(nTotalRays);
    }

    {
        // Gather to master
        List<labelList> allIDs(Pstream::nProcs());
        allIDs[Pstream::myProcNo()] = finishedRayIDs;
        Pstream::gatherList(allIDs);

        List<List<DynamicList<point>>> allPaths(Pstream::nProcs());
        allPaths[Pstream::myProcNo()] = finishedRayPaths;
        Pstream::gatherList(allPaths);

        if (Pstream::master())
        {
            for (label procI = 0; procI < Pstream::nProcs(); ++procI)
            {
                const labelList& ids = allIDs[procI];
                const List<DynamicList<point>>& paths = allPaths[procI];

                forAll(ids, i)
                {
                    const label rayID = ids[i];
                    if (rayID >= 0 && rayID < rayPaths_[laserID].size())
                    {
                        DynamicList<point>& existing =
                            rayPaths_[laserID][rayID];
                        const DynamicList<point>& segment = paths[i];

                        if (existing.empty())
                        {
                            // First segment for this ray
                            existing = segment;
                        }
                        else if (segment.size() > 0)
                        {
                            // Append segment, skipping the first
                            // point if it duplicates the junction
                            label startJ = 0;
                            if
                            (
                                mag(existing.last() - segment[0])
                              < SMALL
                            )
                            {
                                startJ = 1;
                            }
                            for (label j = startJ; j < segment.size(); ++j)
                            {
                                existing.append(segment[j]);
                            }
                        }
                    }
                }
            }
        }
    }

    deposition_.correctBoundaryConditions();

    const scalar TotalQ = fvc::domainIntegrate(deposition_).value();
    Info<< "    Total Q deposited: " << TotalQ << endl;
}


void laserHeatSource::writeRayPathsToVTK()
{
    const Time& runTime = deposition_.time();
    if (Pstream::master()) // Create a directory for the VTK files
    {
        fileName vtkDir;
        if (Pstream::parRun())
        {
            vtkDir = runTime.path()/".."/"VTKs";
        }
        else
        {
            vtkDir = runTime.path()/"VTKs";
        }

        mkDir(vtkDir);

        vtkTimes_.insert(runTime.value());

        forAll(laserNames_, laserI)
        {
            const fileName vtkFileName
            (
                vtkDir/"rays_" + laserNames_[laserI] + "_"
              + Foam::name(runTime.value()) + ".vtk"
            );
            Info<< "Writing " << rayPaths_[laserI].size()
                << " ray paths to " << vtkFileName << endl;
            writeRayPathsToVTK(rayPaths_[laserI], vtkFileName);
        }
    }
}


void laserHeatSource::writeRayPathsToVTK
(
    const List<DynamicList<point>>& rays,
    const fileName& filename
)
{
    OFstream file(filename);

    if (!file.good())
    {
        FatalErrorInFunction
            << "Cannot open file " << filename
            << exit(FatalError);
    }

    label totalPoints = 0;
    label totalLines = 0;

    forAll(rays, rayI)
    {
        totalPoints += rays[rayI].size();
        if (rays[rayI].size() > 1)
        {
            totalLines += (rays[rayI].size() - 1);
        }
    }

    file<< "# vtk DataFile Version 3.0" << nl;
    file<< "Multiple ray data" << nl;
    file<< "ASCII" << nl;
    file<< "DATASET POLYDATA" << nl;

    file<< "POINTS " << totalPoints << " float" << nl;
    forAll(rays, rayI)
    {
        const DynamicList<point>& ray = rays[rayI];
        forAll(ray, pointI)
        {
            const point& pt = ray[pointI];
            file<< pt.x() << " " << pt.y() << " " << pt.z() << nl;
        }
    }

    file<< "LINES " << totalLines << " " << (totalLines * 3) << nl;

    label pointOffset = 0;
    forAll(rays, rayI)
    {
        const DynamicList<point>& ray = rays[rayI];

        for (label i = 0; i < ray.size() - 1; i++)
        {
            file<< "2 " << (pointOffset + i) << " "
                << (pointOffset + i + 1) << nl;
        }

        pointOffset += ray.size();
    }
}


void laserHeatSource::writeRayPathVTKSeriesFile() const
{
    const Time& runTime = deposition_.time();
    if (Pstream::master())
    {
        const fileName vtkDir =
            Pstream::parRun()
          ? runTime.path()/".."/"VTKs"
          : runTime.path()/"VTKs";

        const scalarList times = vtkTimes_.toc();

        SortableList<scalar> sortedTimes(times);

        forAll(laserNames_, laserI)
        {
            const word& laserName = laserNames_[laserI];

            const fileName seriesFile =
                vtkDir/"rays_" + laserName + ".vtk.series";
            Info<< "Writing ray path series file: " << seriesFile << nl;

            OFstream os(seriesFile);
            os.precision(12);

            os  << "{\n"
                << "  \"file-series-version\": \"1.0\",\n"
                << "  \"files\": [\n";

            for (label i = 0; i < sortedTimes.size(); ++i)
            {
                const scalar t = sortedTimes[i];

                os  << "    { \"name\": \"rays_"
                    << laserName << "_" << Foam::name(t)
                    << ".vtk\", \"time\": " << t << " }";

                if (i+1 < sortedTimes.size()) os << ",";
                os << "\n";
            }

            os  << "  ]\n"
                << "}\n";
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
