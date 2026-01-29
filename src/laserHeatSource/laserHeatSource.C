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

#include "laserHeatSource.H"
#include "fvc.H"
#include "constants.H"
#include "findLocalCell.H"
#include "SortableList.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(laserHeatSource, 0);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


void laserHeatSource::createInitialRays
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
) const
{
    DynamicList<vector> initial_points;
    DynamicList<scalar> point_assoc_power;

    const scalarField& yDimI = yDim_;
    const scalar pi = constant::mathematical::pi;

    if (radialPolarHeatSource())
    {
        Info<< "nRadial: " << nRadial << nl
            << "nAngular: "<< nAngular <<endl;

        const scalar rMax = 1.5*beam_radius;
        const label totalSamples = nRadial * nAngular;
        const label samplesPerProc = totalSamples/Pstream::nProcs();
        const label remainder = totalSamples % Pstream::nProcs();
        const label myRank = Pstream::myProcNo();
        const label startIdx = myRank * samplesPerProc + min(myRank, remainder);
        const label endIdx = startIdx + samplesPerProc + (myRank < remainder ? 1 : 0);
        const label localSamples = endIdx - startIdx;

        List<scalar> radialPoints(nRadial);
        for (label iR = 0; iR < nRadial; ++iR)
        {
            // Use sqrt spacing for better Gaussian sampling
            const scalar fraction = scalar(iR + 0.5)/nRadial;
            radialPoints[iR] = rMax * pow(fraction,1.0);
        }

        const point P0 (currentLaserPosition.x(),currentLaserPosition.y(),currentLaserPosition.z());

        // Normalise vector
        const vector V_i(V_incident/(mag(V_incident) + SMALL));

        // // Generate two orthonormal vectors in the plane
        const vector a = (mag(V_i.z()) < 0.9) ? vector(0, 0, 1) : vector(0, 1, 0);
        vector u = (V_i ^ a);
        u = u/mag(u);
        const vector v = (V_i ^ u);
        const vector perturbation (1e-10,1e-10,1e-10);

        for (label localIdx = 0; localIdx < localSamples; ++localIdx)
        {
            const label globalIdx = startIdx + localIdx;

            // Convert global index to angular and radial indices
            const label iTheta = globalIdx/nRadial;
            const label iR = globalIdx % nRadial;

            // Angular discretization
            const scalar theta = 2.0*pi*iTheta/nAngular;

            // Radial discretization (uniform in radius)
            // const scalar r = rMax*(iR + 0.5)/nRadial;
            const scalar r = radialPoints[iR];//if using adaptive sampling

            // Calculate area element
            const scalar deltaTheta = 2.0*pi/nAngular;

            scalar deltaR = 0.0;
            if (iR == 0)
            {
                deltaR = radialPoints[0];
            }
            else
            {
                deltaR = radialPoints[iR] - radialPoints[iR - 1];
            }

            const scalar area = r*deltaR*deltaTheta;

            // Convert to Cartesian coordinates in local plane system
            const scalar x_local = r*cos(theta);
            const scalar y_local = r*sin(theta);

            const vector globalPos = P0 + x_local*u + y_local*v;

            initial_points.append(globalPos + perturbation);

            point_assoc_power.append
            (
                area
               *(
                    Radius_Flavour*Q_cond
                   /(
                        Foam::pow(beam_radius, 2.0)*pi
                    )
                )
               *Foam::exp
                (
                  - Radius_Flavour
                   *(
                        Foam::pow(r, 2.0)/Foam::pow(beam_radius, 2.0)
                    )
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
            // const scalar y_coord = CI[celli].y();
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
                        const point p_1
                        (
                            CI[celli].x()
                          - (yDimI[celli]/2.0)
                          + ((yDimI[celli]/(N_sub_divisions+1))*(Ray_j+1)),
                            CI[celli].y(),
                            CI[celli].z()
                          - (yDimI[celli]/2.0)
                          + ((yDimI[celli]/(N_sub_divisions+1))*(Ray_k+1))
                        );

                        initial_points.append(p_1);

                        point_assoc_power.append
                        (
                            sqr(yDimI[celli]/N_sub_divisions)
                           *(
                                (Radius_Flavour*Q_cond)
                               /(
                                    Foam::pow(beam_radius, 2.0)*pi
                                )
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
    List<pointField> gatheredData(Pstream::nProcs());
    List<scalarField> gatheredData_powers(Pstream::nProcs());

    // Populate and gather the list onto the master processor.
    gatheredData[Pstream::myProcNo()] = initial_points;
    Pstream::gatherList(gatheredData);
    gatheredData_powers[Pstream::myProcNo()] = point_assoc_power;
    Pstream::gatherList(gatheredData_powers);

    // Distibulte the data accross the different processors
    Pstream::broadcastList(gatheredData);
    Pstream::broadcastList(gatheredData_powers);

    // List of initial points
    pointField rayCoords
    (
        ListListOps::combine<Field<vector> >
        (
            gatheredData,
            accessOp<Field<vector> >()
        )
    );

    scalarField rayPowers
    (
        ListListOps::combine<Field<scalar> >
        (
            gatheredData_powers,
            accessOp<Field<scalar> >()
        )
    );

    // Create a list of compactRay objects
    rays.setSize(rayCoords.size());
    forAll(rays, i)
    {
        rays[i] = compactRay
        (
            rayCoords[i], V_incident, rayPowers[i]
        );
        rays[i].globalRayIndex_ = i;
        rays[i].currentCell_ = mesh.findCell(rayCoords[i]);
        rays[i].path_.append(rayCoords[i]);
    }
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
            "Laser_boundary", // rename to laserBoundary?
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
        dimensionedScalar("rayNumber",dimensionSet(0, 0, 0, -0, 0),-1.0)
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
    powderSim_(lookupOrDefault<Switch>("PowderSim", false)),
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
    globalBB_(mesh.bounds())  // Initialize with local bounds first
{
    Info<< "radialPolarHeatSource = " << radialPolarHeatSource_ << endl;

    // Calculate global bounding box
    {
        // Get local bounding box
        boundBox localBB = mesh.bounds();

        // Initialize global bounding box with local bounds
        globalBB_ = localBB;

        // Reduce to get global bounds across all processors
        reduce(globalBB_.min(), minOp<vector>());
        reduce(globalBB_.max(), maxOp<vector>());

        // Inflate the bounding box by 1% to avoid issues with rays starting
        // exactly on the boundary
        globalBB_.inflate(0.01);

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

            laserDicts_.set(laserI, new dictionary(laserEntries[laserI].dict()));

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


    // Update laserBoundary
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

    if (found("HS_lg"))
    {
        FatalErrorInFunction
            << "'HS_lg' is deprecated: please instead specify the laser "
            << "position in time via the laserPositionVsTime sub-dict"
            << exit(FatalError);
    }

    if (found("HS_velocity"))
    {
        FatalErrorInFunction
            << "'HS_velocity' is deprecated: please instead specify the laser "
            << "position in time via the laserPositionVsTime sub-dict"
            << exit(FatalError);
    }

    if (found("HS_Q"))
    {
        FatalErrorInFunction
            << "'HS_Q' is deprecated: please instead specify the laser "
            << "power in time via the laserPowereVsTime sub-dict"
            << exit(FatalError);
    }

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
    // Reset fields
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

            currentLaserPosition[vector::X] += oscAmpX*sin(2*pi*oscFreqX*time);
        }

        // If defined, add oscillation to laser position
        if (dict.found("HS_oscAmpZ"))
        {
            const scalar oscAmpZ(readScalar(dict.lookup("HS_oscAmpZ")));
            const scalar oscFreqZ(readScalar(dict.lookup("HS_oscFreqZ")));
            const scalar pi = constant::mathematical::pi;

            currentLaserPosition[vector::Z] += oscAmpZ*cos(2*pi*oscFreqZ*time);
        }

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

        const label nRadial
        (
            dict.lookupOrDefault<label>("nRadial", 5)
        );

        const label nAngular
        (
            dict.lookupOrDefault<label>("nAngular", 30)
        );

        const label N_sub_divisions
        (
            dict.lookupOrDefault<label>("N_sub_divisions", 1)
        );
        const vector V_incident(dict.lookup("V_incident"));
        const scalar wavelength(readScalar(dict.lookup("wavelength")));
        const scalar e_num_density(readScalar(dict.lookup("e_num_density")));
        const scalar dep_cutoff(dict.lookupOrDefault<scalar>("dep_cutoff", 0.5));

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
    const fvMesh& mesh  = deposition_.mesh();
    const Time& runTime  = mesh.time();
    const scalarField VI = mesh.V();
    const dimensionedScalar time = runTime.time();
    const scalar pi = constant::mathematical::pi;
    const dimensionedScalar a_cond
    (
        "a_cond", dimensionSet(0, 1, 0, 0, 0), laserRadius
    );
    const dimensionedScalar Q_cond
    (
        "Q_cond", dimensionSet(1, 2, -3, 0, 0), currentLaserPower
    );
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

    if (debug)
    {
        Info<< "useLocalSearch: " << useLocalSearch << nl << nl
            << nl << endl;
    }

    // It is assumed that the laser comes in on top y boundary
    const vector normal_interface(0, 1, 0);
    const scalar beam_radius = a_cond.value();

    // Take a references for efficiency and brevity
    const vectorField& nFilteredI = nFiltered;
    const scalarField& alphaFilteredI = alphaFiltered;

    // Create the initial rays
    List<compactRay> rays;
    createInitialRays
    (
        rays,
        mesh,
        currentLaserPosition,
        laserRadius,
        N_sub_divisions,
        nRadial,
        nAngular,
        V_incident,
        Radius_Flavour,
        Q_cond.value(),
        beam_radius
    );

    // remainingGlobalRays will store the rays that have yet to propagate out of
    // the domain or deposit all of their power
    DynamicList<compactRay> remainingGlobalRays(rays);

    // Reset the ray paths list
    if (Pstream::master())
    {
        rayPaths_[laserID].clear();
        rayPaths_[laserID].setSize(rays.size());
    }

    // Calculate the ray power tolerance as a fraction of the max ray power
    // Rays with a power less than this will be ignored
    scalar rayPowerAbsTol = 0;
    {
        scalar maxRayPower = 0.0;
        forAll(remainingGlobalRays, rayI)
        {
            maxRayPower = max(maxRayPower, remainingGlobalRays[rayI].power_);
        }

        rayPowerAbsTol = rayPowerRelTol*maxRayPower;

        Info<< "    Max ray power = " << maxRayPower << nl
            << "    Ray power relative tolerance = " << rayPowerRelTol << nl
            << "    Ray power absolute tolerance = " << rayPowerAbsTol << endl;
    }

    Info<< "    Number of rays: "<< remainingGlobalRays.size() << endl;

    // Propagate the rays through the domain
    while (remainingGlobalRays.size() > 0)
    {
        // Find all rays on the current processor
        // localRays will store the remaining rays on this processor
        DynamicList<compactRay> localRays;
        forAll(remainingGlobalRays, rayI)
        {
            // Take a reference to the current ray
            compactRay& curRay = remainingGlobalRays[rayI];

            // Check if the ray is within the global bound box and its power is
            // greater than a small fraction of the laser power
            if
            (
                globalBB.contains(curRay.position_)
             && curRay.power_ > rayPowerAbsTol
            )
            {
                const label myCellID =
                    findLocalCell
                    (
                        curRay.position_,
                        curRay.currentCell_,
                        mesh,
                        maxLocalSearch,
                        debug
                    );

                if (myCellID != -1)
                {
                    localRays.append(curRay);
                }
            }
        }

        // Propagate the rays through the domain
        forAll(localRays, rayI)
        {
            // Take a reference to the current ray
            compactRay& curRay = localRays[rayI];

            // Find the cell the ray is currently in
            label myCellID =
                findLocalCell
                (
                    curRay.position_,
                    curRay.currentCell_,
                    mesh,
                    maxLocalSearch,
                    debug
                );

            while (myCellID != -1)
            {
                // Calculate the iterator distance as a fraction of the cell size
                const scalar iterator_distance =
                    (0.25/pi)*pow(VI[myCellID], 1.0/3.0);

                // Move the ray by the iterator distance
                curRay.position_ += iterator_distance*curRay.direction_;

                // Find the new cell
                myCellID =
                    findLocalCell
                    (
                        curRay.position_,
                        curRay.currentCell_,
                        mesh,
                        maxLocalSearch,
                        debug
                    );

                // Update the ray's cellID
                curRay.currentCell_ = myCellID;

                if (myCellID == -1)
                {
                    break;
                }

                // Update the rayQ and rayNumber visualisation fields
                rayQ_[myCellID] += curRay.power_;
                rayNumber_[myCellID] = curRay.globalRayIndex_;

                if (curRay.power_ < SMALL)
                {
                    // Update the ray's path
                    curRay.path_.append(curRay.position_);

                    // End of life for the ray
                    break;
                }
                else if
                (
                    mag(nFilteredI[myCellID]) > 0.5
                 && alphaFilteredI[myCellID] >= dep_cutoff
                )
                {
                    // Interface detected
                    // Deposit a fraction of the power and calculate the reflection

                    // Material / dielectric properties
                    const scalar damping_frequency =
                        plasma_frequency*plasma_frequency
                       *constant::electromagnetic::epsilon0.value()
                       *resistivity_in[myCellID];

                    const scalar e_r =
                        1.0
                      - (
                            sqr(plasma_frequency)
                           /(sqr(angular_frequency) + sqr(damping_frequency))
                        );

                    const scalar e_i =
                        (damping_frequency/angular_frequency)
                       *(
                            sqr(plasma_frequency)
                           /(sqr(angular_frequency) + sqr(damping_frequency))
                        );

                    const scalar ref_index =
                        Foam::sqrt
                        (
                            (Foam::sqrt(e_r*e_r + e_i*e_i) + e_r)/2.0
                        );

                    const scalar ext_coefficient =
                        Foam::sqrt
                        (
                            (Foam::sqrt(e_r*e_r + e_i*e_i) - e_r)/2.0
                        );

                    // Incidence angle with consistent normal orientation

                    // Unit surface normal
                    vector n = nFilteredI[myCellID];
                    n /= mag(n);

                    // Incoming direction (towards the interface)
                    vector d = curRay.direction_;
                    const scalar dMag = mag(d);
                    if (dMag <= SMALL)
                    {
                        // Degenerate ray; stop it
                        curRay.power_ = 0.0;

                        // Update the ray's path
                        curRay.path_.append(curRay.position_);

                        // End of life for the ray
                        break;
                    }
                    else
                    {
                        d /= dMag;              // unit propagation direction
                        vector kin = -d;        // direction of incoming wave

                        scalar cosTheta = kin & n;

                        // Flip normal to ensure cosTheta >= 0 (normal into
                        // incident medium)
                        if (cosTheta < 0.0)
                        {
                            n = -n;
                            cosTheta = -cosTheta;
                        }

                        // Clamp to [0,1] to avoid NaN from acos
                        cosTheta =
                            Foam::max
                            (
                                Foam::min(cosTheta, scalar(1.0)),
                                scalar(0.0)
                            );

                        const scalar theta_in = std::acos(cosTheta);

                        // Fresnel-like optics

                        const scalar sinTheta = Foam::sin(theta_in);

                        const scalar alpha_laser =
                            Foam::sqrt
                            (
                                Foam::sqrt
                                (
                                    sqr
                                    (
                                        sqr(ref_index)
                                      - sqr(ext_coefficient)
                                      - sqr(sinTheta)
                                    )
                                  + 4.0*sqr(ref_index)*sqr(ext_coefficient)
                                )
                              + sqr(ref_index)
                              - sqr(ext_coefficient)
                              - sqr(sinTheta)/2.0
                            );

                        const scalar beta_laser =
                            Foam::sqrt
                            (
                                (
                                    Foam::sqrt
                                    (
                                        sqr
                                        (
                                            sqr(ref_index)
                                          - sqr(ext_coefficient)
                                          - sqr(sinTheta)
                                        )
                                      + 4.0*sqr(ref_index)*sqr(ext_coefficient)
                                    )
                                  - sqr(ref_index)
                                  + sqr(ext_coefficient)
                                  + sqr(sinTheta)
                                )/2.0
                            );

                        const scalar cosTheta_in = Foam::cos(theta_in);

                        scalar R_s =
                            (
                                sqr(alpha_laser)
                              + sqr(beta_laser)
                              - 2.0*alpha_laser*cosTheta_in
                              + sqr(cosTheta_in)
                            )
                           /(
                                sqr(alpha_laser)
                              + sqr(beta_laser)
                              + 2.0*alpha_laser*cosTheta_in
                              + sqr(cosTheta_in)
                            );

                        scalar R_p =
                            R_s
                           *(
                                sqr(alpha_laser)
                              + sqr(beta_laser)
                              - 2.0*alpha_laser*sinTheta*Foam::tan(theta_in)
                              + sqr(sinTheta)*sqr(Foam::tan(theta_in))
                            )
                           /(
                                sqr(alpha_laser)
                              + sqr(beta_laser)
                              + 2.0*alpha_laser*sinTheta*Foam::tan(theta_in)
                              + sqr(sinTheta)*sqr(Foam::tan(theta_in))
                            );

                        // Clamp reflectivities and absorptivity to [0,1]

                        R_s =
                            Foam::max(Foam::min(R_s, scalar(1.0)), scalar(0.0));
                        R_p =
                            Foam::max(Foam::min(R_p, scalar(1.0)), scalar(0.0));

                        scalar R = 0.5*(R_s + R_p);
                        scalar absorptivity = 1.0 - R;

                        absorptivity =
                            Foam::max
                            (
                                Foam::min(absorptivity, scalar(1.0)),
                                scalar(0.0)
                            );

                        if (debug)
                        {
                            Info<< "ray " << curRay.globalRayIndex_
                                << ", cell " << myCellID
                                << ", theta_in = " << theta_in
                                << ", R_s = " << R_s
                                << ", R_p = " << R_p
                                << ", absorptivity = " << absorptivity << endl;
                        }

                        // Deposit and reflect

                        deposition_[myCellID] +=
                            absorptivity*curRay.power_/VI[myCellID];

                        curRay.power_ *= (1.0 - absorptivity);

                        // Specular reflection
                        // reflect original direction about n
                        // n points into incident medium
                        const vector dRef =
                            curRay.direction_ - 2.0*(curRay.direction_ & n)*n;

                        curRay.direction_ = dRef;
                    }
                }
                else if (alphaFilteredI[myCellID] >= dep_cutoff)
                {
                    // Bulk metal
                    // Assume fully absorbing, no further propagation

                    if (debug)
                    {
                        Info<< "Bulk absorption at cell " << myCellID
                            << ", alpha = " << alphaFilteredI[myCellID]
                            << ", |n| = " << mag(nFilteredI[myCellID])
                            << ", power = " << curRay.power_ << endl;
                    }

                    // Deposit all remaining ray power in this cell
                    deposition_[myCellID] += curRay.power_/VI[myCellID];

                    // Kill the ray
                    curRay.power_ = 0.0;

                    // Update the ray's path
                    curRay.path_.append(curRay.position_);

                    // End of life for the ray
                    break;
                }

                // Update the ray's path
                curRay.path_.append(curRay.position_);
            }
        }

        // Sync all remaining local rays globally so remainingGlobalRays will
        // be the same on all processors
        remainingGlobalRays = localRays;
        Pstream::combineGather(remainingGlobalRays, combineRayLists());
        Pstream::broadcast(remainingGlobalRays);

        // Record the latest ray paths
        // Note that once a ray has left the domain then its global path is no
        // longer updated so its path will be the final full path
        if (Pstream::master())
        {
            forAll(remainingGlobalRays, rI)
            {
                const compactRay& curRay = remainingGlobalRays[rI];
                const label rayID = curRay.globalRayIndex_;
                rayPaths_[laserID][rayID] = curRay.path_;
            }
        }
    }

    deposition_.correctBoundaryConditions();

    const scalar TotalQ = fvc::domainIntegrate(deposition_).value();
    Info<< "    Total Q deposited: " << TotalQ <<endl;
}


void laserHeatSource::writeRayPathsToVTK()
{
    const Time& runTime = deposition_.time();
    if (Pstream::master())
    {
        // Create a directory for the VTK files
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

        // Record this time as having a VTK file
        vtkTimes_.insert(runTime.value());

        // Write ray paths to VTK files
        forAll(laserNames_, laserI)
        {
            const fileName vtkFileName
            (
                vtkDir/"rays_" + laserNames_[laserI] + "_"
              + Foam::name(runTime.value()) + ".vtk"
            );
            Info<< "Writing " << rayPaths_[laserI].size() << " ray paths to "
                << vtkFileName << endl;
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

    // Calculate total number of points and lines
    label totalPoints = 0;
    label totalLines = 0;

    forAll(rays, rayI)
    {
        totalPoints += rays[rayI].size();
        if (rays[rayI].size() > 1)
        {
            // n-1 line segments per ray
            totalLines += (rays[rayI].size() - 1);
        }
    }

    // VTK header
    file<< "# vtk DataFile Version 3.0" << nl;
    file<< "Multiple ray data" << nl;
    file<< "ASCII" << nl;
    file<< "DATASET POLYDATA" << nl;

    // Write all points
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

    // Write line connectivity
    // Each line segment is defined by 2 points
    file<< "LINES " << totalLines << " " << (totalLines * 3) << nl;

    label pointOffset = 0;
    forAll(rays, rayI)
    {
        const DynamicList<point>& ray = rays[rayI];

        // Connect consecutive points within this ray only
        for (label i = 0; i < ray.size() - 1; i++)
        {
            file<< "2 " << (pointOffset + i) << " " << (pointOffset + i + 1)
                << nl;
        }

        // Update offset for next ray
        pointOffset += ray.size();
    }
}


void laserHeatSource::writeRayPathVTKSeriesFile() const
{
    const Time& runTime = deposition_.time();
    if (Pstream::master())
    {
        // Directory that already holds the legacy .vtk files
        const fileName vtkDir =
            Pstream::parRun()
          ? runTime.path()/".."/"VTKs"
          : runTime.path()/"VTKs";

        // Get the list of times when VTK files were written
        const scalarList times = vtkTimes_.toc();

        // Ensure chronological order
        SortableList<scalar> sortedTimes(times);

        // Write ray paths to VTK files
        forAll(laserNames_, laserI)
        {
            const word& laserName = laserNames_[laserI];

            const fileName seriesFile = vtkDir/"rays_" + laserName + ".vtk.series";
            Info<< "Writing ray path series file: " << seriesFile << nl;

            OFstream os(seriesFile);
            os.precision(12); // good numeric precision for times

            os  << "{\n"
                << "  \"file-series-version\": \"1.0\",\n"
                << "  \"files\": [\n";

            for (label i = 0; i < sortedTimes.size(); ++i)
            {
                const scalar t = sortedTimes[i];

                // Filenames follow your convention:
                //   rays_<laserName>_<timeName>.vtk
                // We reconstruct <timeName> with Foam::name(t).
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
