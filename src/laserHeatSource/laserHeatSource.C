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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(laserHeatSource, 0);


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
    refineFlag_
    (
        IOobject
        (
            "refineflag", // rename refineFlag
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("refineflag", dimensionSet(0,0,0,0,0), 0.0)
    ),
    powderSim_(lookupOrDefault<Switch>("PowderSim", false)),
    Radial_Polar_HS_(lookupOrDefault<Switch>("Radial_Polar_HS", false)),
    laserNames_(0),
    laserDicts_(0),
    timeVsLaserPosition_(0),
    timeVsLaserPower_(0)
{
    // Initialise the laser power and position
    if (found("lasers"))
    {
        const PtrList<entry> laserEntries(lookup("lasers"));

        laserNames_.setSize(laserEntries.size());
        laserDicts_.setSize(laserEntries.size());
        timeVsLaserPosition_.setSize(laserEntries.size());
        timeVsLaserPower_.setSize(laserEntries.size());

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
bool laserHeatSource::pointOnFace(const face& f, const point& p, const pointField& meshPoints){

    pointField fPoints(f.points(meshPoints));

    // Step 1: Check if point lies on the plane of the face
    vector normal = f.normal(meshPoints);
    scalar d = - (normal & fPoints[0]); // Plane equation: n·x + d = 0
    scalar distToPlane = (normal & p) + d;
if(debug){
     Info<<"dist to plane: "<<distToPlane<<endl;
}
    if (mag(distToPlane) > SMALL){  // Point is not in the plane
        return false;
        }
        else{
            return true;
        }




}

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

        Info<< "Laser " << laserNames_[laserI] << nl
            << "Laser mean position = " << currentLaserPosition << nl
            << "Laser power = " << currentLaserPower << endl;

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

        Info<< "Laser position including any oscillation = "
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

        updateDeposition
        (
            alphaFiltered,
            nFiltered,
            resistivity_in,
            laserNames_[laserI],
            currentLaserPosition,
            currentLaserPower,
            laserRadius,
            N_sub_divisions,
            V_incident,
            wavelength,
            e_num_density,
            dep_cutoff,
            Radius_Flavour,
            useLocalSearch,
            maxLocalSearch
        );
    }
}

void laserHeatSource::updateDeposition
(
    const volScalarField& alphaFiltered,
    const volVectorField& nFiltered,
    const volScalarField& resistivity_in,
    const word& laserName,
    const vector& currentLaserPosition,
    const scalar currentLaserPower,
    const scalar laserRadius,
    const label N_sub_divisions,
    const vector& V_incident,
    const scalar wavelength,
    const scalar e_num_density,
    const scalar dep_cutoff,
    const scalar Radius_Flavour,
    const Switch useLocalSearch,
    const label maxLocalSearch
)
{
    const fvMesh& mesh  = deposition_.mesh();
    const Time& runTime  = mesh.time();
    const dimensionedScalar time = runTime.time();
    const dimensionedScalar pi = constant::mathematical::pi;
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
        2.0*pi.value()*constant::universal::c.value()/wavelength;

    if (debug)
    {
        Info<< "useLocalSearch: " << useLocalSearch << nl << nl
            << nl << endl;
    }

    // It is assumed that the laser comes in on top y boundary
    const vector normal_interface(0, 1, 0);

    const scalar beam_radius =
        a_cond.value();
    //    /Foam::cos
    //     (
    //         Foam::acos
    //         (
    //             (normal_interface & (V_incident/mag(V_incident)))
    //            /(mag(normal_interface)*mag(V_incident/mag(V_incident)))
    //         )
    //     );

    // Adjust sample radius for if beam is not normal too top boundary

    const scalar CosTheta_incident =
        Foam::cos
        (
            Foam::acos
            (
                (normal_interface & (V_incident/mag(V_incident)))
               /(mag(normal_interface)*mag(V_incident/mag(V_incident)))
            )
        );

    // if (debug)
    // {
    //     Info<< "cos(theta): " << CosTheta_incident << endl;
    // }

    scalar listLength(0);
    DynamicList<vector> initial_points(listLength, vector::zero);
    initial_points.clear();

    // List with size equal to number of processors
    List<pointField> gatheredData1(Pstream::nProcs());

    // Take a references for efficiency and brevity
    const vectorField& CI = mesh.C();
    const scalarField& yDimI = yDim_;
    const vectorField& nFilteredI = nFiltered;
    const scalarField& alphaFilteredI = alphaFiltered;









 





    
    if(Radial_Polar_HS()==true){

        label nRings = 200;
        label nAngles = 120;
        scalar rMax = 1.5*beam_radius;
        point P0 (currentLaserPosition.x(),currentLaserPosition.y(),currentLaserPosition.z());
        // Info<<" \n \n TESTTTTTTTT \n \n"<<endl;
        vector V_i(V_incident/mag(V_incident)); //normalise vector in-case user hasnt

        // Generate two orthonormal vectors in the plane
        vector a = (mag(V_i.z()) < 0.9) ? vector(0, 0, 1) : vector(0, 1, 0);
        vector u = (V_i ^ a);
        u = u/mag(u);
        vector v = (V_i ^ u);
        vector perturbation (1e-9,1e-9,1e-9);

        // if(mesh.findCell(P0)!=-1){
        initial_points.append(P0 + perturbation);
        // }

        for (label i = 1; i < nRings; ++i)
        {
        scalar r = rMax * scalar(i) / scalar(nRings);  // linear spacing
        for (label j = 0; j < nAngles; ++j)

            {

            scalar theta = 2.0 * pi.value() * scalar(j) / scalar(nAngles);

            vector offset = r * (cos(theta) * u + sin(theta) * v);

            
            // Info<<"test"<<mesh.findCell(P0 + offset)<<endl;
            // if(mesh.findCell(P0 + offset)!=-1){
            
    //         pointslistGlobal1[i] += perturbation;
                initial_points.append(P0 + offset + perturbation);
            // }
            // Info<<"i disc: "<<i<<", j disc: "<<j<<endl;
            // Info<<P0 + offset<<endl;
            }

        }

        // Info<<initial_points<<endl;
        // Info<<initial_points[270]<<endl;


    //     forAll(initial_points, i)
    // {
    //     //find if a point is on a face and add a small perturbation
    //     const pointField& points = mesh.points();

    //     forAll(mesh.faces(),faceI){
    //             // Info<<"FACE: "<<faceI<<endl;
    //             const face& f = mesh.faces()[faceI];

    //                     if (pointOnFace(f, initial_points[i], points))
    //     {
    //         if(debug){
    //         Info << "Point is on face " << faceI << nl;
    //         }
    //         vector perturbation (1e-9,1e-9,1e-9);
    //         initial_points[i] += perturbation;

    //     }

                
    //     }
   
    // }

    
    }
    else{

           forAll(CI, celli)
    {
        const scalar x_coord = CI[celli].x();
        // const scalar y_coord = CI[celli].y();
        const scalar z_coord = CI[celli].z();

        if
        (
            (
                Foam::pow(x_coord - currentLaserPosition.x(), 2.0)
              + Foam::pow(z_coord - currentLaserPosition.z(), 2.0)
             <= Foam::pow(1.5*beam_radius, 2.0)
            )
         && (laserBoundary_[celli] > SMALL)
        )
        {
            // rayNumber_[celli] = 1.0;
            refineFlag_[celli] += 0.5;
            // initial_points.append(CI[celli]);

            for (label Ray_j = 0; Ray_j < N_sub_divisions; Ray_j++)
            {
                for(label Ray_k = 0; Ray_k < N_sub_divisions; Ray_k++)
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
                    initial_points.append(p_1);
                }
            }

        }
    }
    }










    // Populate and gather the list onto the master processor.
    gatheredData1[Pstream::myProcNo()] = initial_points;
    Pstream::gatherList(gatheredData1);

    // Distibulte the data accross the different processors
    Pstream::scatterList(gatheredData1);

    // List of initial points
    pointField pointslistGlobal1
    (
        ListListOps::combine<Field<vector> >
        (
            gatheredData1,
            accessOp<Field<vector> >()
        )
    );
    // Info<<"synched points: "<<pointslistGlobal1<<endl;






///////////////////////////////////////////////////need to loop over all global points and add perturbation if they are exactly on a face 

    //         forAll(pointslistGlobal1, i)
    // {
    //     //find if a point is on a face and add a small perturbation
    //     const pointField& points = mesh.points();

    //     forAll(mesh.faces(),faceI){
    //             // Info<<"FACE: "<<faceI<<endl;
    //             const face& f = mesh.faces()[faceI];

    //                     if (pointOnFace(f, pointslistGlobal1[i], points))
    //     {
    //         if(debug){
    //         Info << "Point is on face " << faceI << nl;
    //         }
    //         vector perturbation (1e-9,1e-9,1e-9);
    //         pointslistGlobal1[i] += perturbation;

    //     }

                
    //     }
   
    // }
///////////////////////////////////////////////////need to loop over all global points and add perturbation if they are exactly on a face 



    
    // For each beam, store the starting point and locations at which the rays
    // change direction. Also, store the global ordered index of the ray
    // direction-change points
    PtrList<DynamicList<vector>> beamDirectionChangePoints
    (
        pointslistGlobal1.size()
    );
    PtrList<DynamicList<int>> beamDirectionChangeOrder
    (
        pointslistGlobal1.size()
    );

    // Initialise: beams will likely change direction less than 100 times
    forAll(beamDirectionChangePoints, rayI)
    {
        beamDirectionChangePoints.set
        (
            rayI,
            new DynamicList<vector>(100)
        );
        beamDirectionChangeOrder.set
        (
            rayI,
            new DynamicList<int>(100)
        );

        // Add initial point
        beamDirectionChangePoints[rayI].append(pointslistGlobal1[rayI]);
        beamDirectionChangeOrder[rayI].append(0);
    }

    // Store the list of cell indices where the ray tips are located; these will
    // be used by the the findLocalSearch function when looking for the new tip
    // cell indices
    labelList rayCellIDs(pointslistGlobal1.size(), -1);

    const point DUMMYMAX(-GREAT, -GREAT, -GREAT);
    const scalar DUMMYSCAL(-GREAT);

    // Loop over all starting points
    Info<<"Calculating laser beam rays" << endl;
    forAll(pointslistGlobal1, i)
    {
        if (debug)
        {
            Info<< "Beam " << i << endl;
        }

        vector V2(V_incident/mag(V_incident));
        point V1_tip(pointslistGlobal1[i]);

        const point mid
        (
            currentLaserPosition.x(),
            pointslistGlobal1[i].y(),
            currentLaserPosition.z()
        );

        if(Radial_Polar_HS()==true){
            const point mid
            (
            currentLaserPosition.x(),
            currentLaserPosition.y(),//pointslistGlobal1[i].y(),
            currentLaserPosition.z()
            );
        }
        // else{

        // }

        const vector x1 = mid - 10.0*V2;
        const vector x2 = mid + 10.0*V2;
        const vector x0
        (
            pointslistGlobal1[i].x(),
            pointslistGlobal1[i].y(),
            pointslistGlobal1[i].z()
        );
        // Info<<"HERE"<<endl;
        // Cross product to find distance to beam central axis
        const scalar dist = mag(((x0 - x1)^(x0 - x2)))/mag(x2 - x1);

        // Global index to track the order of the ray direction-changes
        // This is only used for post-processing to write VTKs of the beams
        label directionChangeOrderI = 0;

        // Info<<"HERE2"<<endl;
        

        scalar Q =
            (
                CosTheta_incident/(N_sub_divisions*N_sub_divisions)
            )
           *
           (
               (Radius_Flavour*Q_cond.value())
              /(
                  Foam::pow(a_cond.value(), 2.0)*pi.value()
               )
           )
          *Foam::exp
           (
             - Radius_Flavour
              *(
                  Foam::pow(dist, 2.0)/Foam::pow(a_cond.value(), 2.0)
               )
           );

if(Radial_Polar_HS()==true){
                Q =
           (
               (Radius_Flavour*Q_cond.value())
              /(
                  pointslistGlobal1.size()*Foam::pow(a_cond.value(), 2.0)*pi.value()
               )
           )
          *Foam::exp
           (
             - Radius_Flavour
              *(
                  Foam::pow(dist, 2.0)/Foam::pow(a_cond.value(), 2.0)
               )
           );
       
}

            // Info<<"dist: "<<dist<<", Q: "<<Q<<endl;
           


        // ID of the processor that contains the beam tip
        label tipProcID = -1;

        while (Q > 1.0e-9)//getting stuck in this loop - newHS
        {
            // Info<<"HERE4"<<Q<<endl;
            // Track when the tip changes direction for post-processing the rays
            bool beamChangedDirection = false;

            // Search for the cell that contains the local beam tip
            // Only the processor that contained the old tip will perform the
            // search, or all processor will search if the old tip is not on any
            // processor
            label myCellId = -1;
            if (tipProcID == Pstream::myProcNo() || tipProcID == -1)
            {
                if (useLocalSearch)
                {
                    myCellId =
                        findLocalCell
                        (
                            V1_tip, rayCellIDs[i], mesh, maxLocalSearch, debug
                        );
                }
                else
                {
                    myCellId = mesh.findCell(V1_tip);
                }
            }
            // Info<<myCellId<<"\t"<<Q<<endl;
            // Q/=mesh.V()[myCellId];
            // Proc ID where the tip is located
            // If the tip in not on any processor, then this is set to -1
            if (myCellId != -1)
            {
                tipProcID = Pstream::myProcNo();
            }
            else
            {
                tipProcID = -1;
            }
            reduce(tipProcID, maxOp<label>());
            // Info<<"HERE 1: "<<myCellId<<endl;
            if (myCellId != -1)
            {
                // Set test field to beam flavour
                rayNumber_[myCellId] = i + 1;
                rayQ_[myCellId] += Q;
                // Info<<"HERE 2: "<<myCellId<<endl;
                if
                (
                    mag(nFilteredI[myCellId]) > 0.5
                 && alphaFilteredI[myCellId] >= dep_cutoff
                )
                {
                    const scalar damping_frequency =
                        plasma_frequency*plasma_frequency
                       *constant::electromagnetic::epsilon0.value()
                       *resistivity_in[myCellId];

                    const scalar e_r =
                        1.0
                      - (
                            sqr(plasma_frequency)/(sqr(angular_frequency)
                          + sqr(damping_frequency))
                        );

                    const scalar e_i =
                        (damping_frequency/angular_frequency)
                        *(
                            plasma_frequency*plasma_frequency
                           /(
                                angular_frequency*angular_frequency
                              + damping_frequency*damping_frequency
                            )
                        );

                    const scalar ref_index =
                        Foam::sqrt
                        (
                            (Foam::sqrt((e_r*e_r) +(e_i*e_i)) + e_r)/2.0
                        );

                    const scalar ext_coefficient =
                        Foam::sqrt
                        (
                            (Foam::sqrt((e_r*e_r) +(e_i*e_i)) - e_r)/2.0
                        );

                    scalar argument =
                        (
                            V2 & nFilteredI[myCellId]
                        )/(mag(V2)*mag(nFilteredI[myCellId]));
                    // Info<<"HERE 2"<<endl;
                    if (argument >= (1.0 - SMALL))
                    {
                        argument = 1.0;
                    }
                    else if (argument <= (-1.0 + SMALL))
                    {
                        argument = -1.0;
                    }

                    const scalar theta_in = std::acos(argument);

                    const scalar alpha_laser =
                        Foam::sqrt
                        (
                            Foam::sqrt
                            (
                                sqr
                                (
                                    sqr(ref_index)
                                  - sqr(ext_coefficient)
                                  - sqr(Foam::sin(theta_in))
                                )
                              + (
                                    4.0*sqr(ref_index)*sqr(ext_coefficient)
                                )
                            )
                          + sqr(ref_index)
                          - sqr(ext_coefficient)
                          - sqr(Foam::sin(theta_in))/2.0
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
                                      - sqr(Foam::sin(theta_in))
                                    )
                                  + 4.0*sqr(ref_index)*sqr(ext_coefficient)
                                )
                              - sqr(ref_index)
                              + sqr(ext_coefficient)
                              + sqr(Foam::sin(theta_in))
                            )/2.0
                        );

                    const scalar R_s =
                        (
                            (
                                sqr(alpha_laser)
                              + sqr(beta_laser)
                              - 2.0*alpha_laser*Foam::cos(theta_in)
                              + sqr(Foam::cos(theta_in))
                            )
                           /(
                               sqr(alpha_laser)
                             + sqr(beta_laser)
                             + 2.0*alpha_laser*Foam::cos(theta_in)
                             + sqr(Foam::cos(theta_in))
                           )
                        );

                    const scalar R_p =
                        R_s
                       *(
                           (
                               sqr(alpha_laser)
                             + sqr(beta_laser)
                             - (
                                   2.0*alpha_laser*Foam::sin(theta_in)
                                  *Foam::tan(theta_in)
                               )
                             + (
                                   sqr(Foam::sin(theta_in))
                                  *sqr(Foam::tan(theta_in))
                               )
                           )
                          /(
                              sqr(alpha_laser)
                            + sqr(beta_laser)
                            + (
                                  2.0*alpha_laser*Foam::sin(theta_in)
                                 *Foam::tan(theta_in)
                              )
                            + sqr(Foam::sin(theta_in))*sqr(Foam::tan(theta_in))
                           )
                       );
                    const scalar absorptivity = 1.0 - ((R_s + R_p)/2.0);

                    // Info<<i<<"\t"<<absorptivity<<endl;

                    // Sometimes the ray can be reflected and 'skip' along the
                    // interface cells - this is unphysical and the ray should
                    // traverse  without depositing any energy so set Q to 0 in
                    // this instance
                    if (theta_in >= pi.value()/2.0)
                    {
                        Q *= 0.0;
                        deposition_[myCellId] += absorptivity*Q/yDimI[myCellId];
                        if (debug)
                        {
                            errorTrack_[myCellId] -= 1.0;
                        }
                        beamChangedDirection = true;
                    }
                    // else{}
                    else
                    {
                        deposition_[myCellId] += absorptivity*Q/yDimI[myCellId];
                        Q *= 1.0 - absorptivity;
                        V2 -=
                            (
                                (
                                    (
                                        ((2.0*V2) & nFilteredI[myCellId])
                                       /(
                                            mag(nFilteredI[myCellId])
                                           *mag(nFilteredI[myCellId])
                                        )
                                    )
                                )*nFilteredI[myCellId]
                            );

                        beamChangedDirection = true;
                    }
                }
                else
                {
                    // If the ray step size happens to be large enough that it
                    // skips through the interface send ray back the way it came
                    if
                    (
                        alphaFilteredI[myCellId] > dep_cutoff
                     && mag(nFilteredI[myCellId]) < 0.5
                    )
                    {
                        const scalar damping_frequency =
                            plasma_frequency*plasma_frequency
                           *constant::electromagnetic::epsilon0.value()
                           *resistivity_in[myCellId];

                        const scalar e_r =
                            1.0
                          - (
                                sqr(plasma_frequency)/(sqr(angular_frequency)
                              + sqr(damping_frequency))
                            );

                        const scalar e_i =
                            (damping_frequency/angular_frequency)
                            *(
                                 plasma_frequency*plasma_frequency
                                /(
                                     angular_frequency*angular_frequency
                                   + damping_frequency*damping_frequency
                                 )
                            );

                        const scalar ref_index =
                            Foam::sqrt
                            (
                                (Foam::sqrt((e_r*e_r) +(e_i*e_i)) + e_r)/2.0
                            );

                        const scalar ext_coefficient =
                            Foam::sqrt
                            (
                                (Foam::sqrt((e_r*e_r) +(e_i*e_i)) - e_r)/2.0
                            );

                        if (debug)
                        {
                            errorTrack_[myCellId] += 1.0;
                        }

                        const scalar theta_in = 0.0;

                        const scalar alpha_laser =
                            Foam::sqrt
                            (
                                (
                                    Foam::sqrt
                                    (
                                        sqr
                                        (
                                            sqr(ref_index)
                                          - sqr(ext_coefficient)
                                          - sqr(Foam::sin(theta_in))
                                        )
                                      + 4.0*sqr(ref_index)*sqr(ext_coefficient)
                                    )
                                  + sqr(ref_index)
                                  - sqr(ext_coefficient)
                                  - sqr(Foam::sin(theta_in))
                                )/2.0
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
                                          - sqr(Foam::sin(theta_in))
                                        )
                                      + 4.0*sqr(ref_index)*sqr(ext_coefficient)
                                    )
                                  - sqr(ref_index)
                                  + sqr(ext_coefficient)
                                  + sqr(Foam::sin(theta_in))
                                )/2.0
                            );

                        const scalar R_s =
                            (
                                (
                                    sqr(alpha_laser)
                                  + sqr(beta_laser)
                                  - (2.0*alpha_laser*Foam::cos(theta_in))
                                  + sqr(Foam::cos(theta_in))
                                )
                               /(
                                   sqr(alpha_laser)
                                 + sqr(beta_laser)
                                 + (2.0*alpha_laser*Foam::cos(theta_in))
                                 + sqr(Foam::cos(theta_in))
                                )
                            );

                        const scalar R_p =
                            R_s
                           *(
                               (
                                   sqr(alpha_laser)
                                 + sqr(beta_laser)
                                 - 2.0*alpha_laser*Foam::sin(theta_in)
                                  *Foam::tan(theta_in)
                                 + sqr(Foam::sin(theta_in))
                                  *sqr(Foam::tan(theta_in))
                               )
                              /(
                                  sqr(alpha_laser)
                                + sqr(beta_laser)
                                + 2.0*alpha_laser*Foam::sin(theta_in)
                                 *Foam::tan(theta_in)
                                + sqr(Foam::sin(theta_in))
                                 *sqr(Foam::tan(theta_in))
                               )
                           );

                        const scalar absorptivity = 1.0 - ((R_s + R_p)/2.0);

                        // If the ray slips through the interface (unlikely)
                        // send it back the way it came because it must have
                        // been at 0 degrees anyway
                        V2 = -V2;

                        beamChangedDirection = true;
                        deposition_[myCellId] += absorptivity*Q/yDimI[myCellId];
                        Q *= (1.0 - absorptivity);
                    }
                }
             }
             else
             {
                 // The tip is not on this processor for one of two reasons:
                 // 1. the tip left the entire global domain
                 // 2. the tip is on another processor
                 V2 = DUMMYMAX;
                 Q = DUMMYSCAL;
                 beamChangedDirection = true;
             }
// Info<<"HERE5"<<endl;
             reduce(V2, maxOp<vector>());
             reduce(Q, maxOp<scalar>());

             // Update seed cells for local search
             rayCellIDs[i] = myCellId;
// Info<<"HERE56"<<endl;
             if (tipProcID == Pstream::myProcNo())
             {
                 label myCellIdnext =
                     findLocalCell
                     (
                         V1_tip, rayCellIDs[i], mesh, maxLocalSearch, debug
                     );

                 if (myCellIdnext != -1)
                 {
                     while (myCellIdnext == myCellId)
                     {
                         if (beamChangedDirection)
                         {
                             // Write current tip position to array
                             beamDirectionChangePoints[i].append(V1_tip);
                             beamDirectionChangeOrder[i].append
                             (
                                 directionChangeOrderI
                             );
                             beamChangedDirection = false;
                         }

                         const scalar iterator_distance =
                             (0.5/pi.value())*yDimI[myCellId];

                         if (debug)
                         {
                             Info<<"iterator_distance = "
                                 << iterator_distance << endl;
                         }

                         V1_tip += iterator_distance*V2;
                         myCellIdnext =
                             findLocalCell
                             (
                                 V1_tip,
                                 rayCellIDs[i],
                                 mesh,
                                 maxLocalSearch,
                                 debug
                             );
                     }
                 }
                 else
                 {
                     V1_tip = DUMMYMAX;
                 }
                // Info<<"HERE567"<<endl;
                 // Update direction-change ordered index
                 directionChangeOrderI++;

                 // Update seed cells for local search
                 rayCellIDs[i] = myCellIdnext;
             }
             else
             {
                 V1_tip = DUMMYMAX;
             }
             reduce(V1_tip, maxOp<vector>());

             if (rayCellIDs[i] == -1)
             {
                 tipProcID = -1;
             }
             reduce(tipProcID, maxOp<label>());
// Info<<"HERE5678: "<<Q<<endl;
             // Sync direction-change ordered index
             reduce(directionChangeOrderI, maxOp<int>());
         };
        //  Info<<"HERE56789"<<endl;
     }


     const scalar TotalQ = fvc::domainIntegrate(deposition_).value();
     Info<< "Total Q deposited this timestep: " << TotalQ <<endl;

     // Combine rays across procs
     if (runTime.outputTime() && Pstream::parRun())
     {
         if (debug)
         {
             Info<< "Parallel syncing beams!" << endl;
         }

         // The ray starting points were added to the beamDirectionChangePoints list
         // on all procs, so we will remove them from all procs apart from the master
         // Note: the beamDirectionChangePoints list is only synced at output times
         // and will only be correct on the master proc which writes them
         if (!Pstream::master())
         {
             forAll(beamDirectionChangePoints, rayI)
             {
                 beamDirectionChangePoints[rayI] =
                     SubField<vector>
                     (
                         beamDirectionChangePoints[rayI],
                         beamDirectionChangePoints[rayI].size() - 1,
                         1
                     );

                 beamDirectionChangeOrder[rayI] =
                     SubField<int>
                     (
                         beamDirectionChangeOrder[rayI],
                         beamDirectionChangeOrder[rayI].size() - 1,
                         1
                     );
             }
         }

         // Sync beams across procs
         forAll(beamDirectionChangePoints, rayI)
         {
             {
                 List<List<vector>> gatheredField(Pstream::nProcs());
                 gatheredField[Pstream::myProcNo()] =
                     beamDirectionChangePoints[rayI];
                 Pstream::gatherList(gatheredField);

                 beamDirectionChangePoints[rayI] =
                     ListListOps::combine<List<vector>>
                     (
                         gatheredField,
                         accessOp<List<vector>>()
                     );
             }

             {
                 List<List<int>> gatheredField(Pstream::nProcs());
                 gatheredField[Pstream::myProcNo()] =
                     beamDirectionChangeOrder[rayI];
                 Pstream::gatherList(gatheredField);

                 beamDirectionChangeOrder[rayI] =
                     ListListOps::combine<List<int>>
                     (
                         gatheredField,
                         accessOp<List<int>>()
                     );
             }

             // Re-order the list
             if (Pstream::master())
             {
                 SortableList<int> sortedOrder(beamDirectionChangeOrder[rayI]);
                 List<vector> unsortedPoints(beamDirectionChangePoints[rayI]);
                 List<int> unsortedOrder(beamDirectionChangeOrder[rayI]);
                 forAll(sortedOrder, i)
                 {
                     beamDirectionChangePoints[rayI][i] =
                         unsortedPoints[sortedOrder.indices()[i]];
                     beamDirectionChangeOrder[rayI][i] =
                         unsortedOrder[sortedOrder.indices()[i]];
                 }
             }
         }
     }


     // Write rays
     if (runTime.outputTime() && Pstream::master())
     {
         if (debug)
         {
             forAll(beamDirectionChangePoints, rayI)
             {
                 Info<< "ray " << rayI << endl;
                 forAll(beamDirectionChangePoints[rayI], i)
                 {
                     Info<< "    " << beamDirectionChangePoints[rayI][i]
                         << endl;
                 }
             }
         }

         // Write rays in VTK format
         // See
         // https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html

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

         // Create a VTK file
         OFstream rayVtkFile
         (
             vtkDir/laserName + "_rays_"
           + Foam::name(runTime.timeIndex()) + ".vtk"
         );

         Info<< "Writing rays to " << rayVtkFile.name() << endl;

         // Write header
         rayVtkFile
             << "# vtk DataFile Version 2.0" << nl
             << "Rays" << nl
             << "ASCII" << endl;

         // Count the number of points and calculate the offset for each ray
         label nRayPoints = 0;
         labelList pointIdOffset(beamDirectionChangePoints.size(), 0);
         forAll(beamDirectionChangePoints, rayI)
         {
             nRayPoints += beamDirectionChangePoints[rayI].size();

             if (rayI > 0)
             {
                 pointIdOffset[rayI] =
                     pointIdOffset[rayI - 1]
                   + beamDirectionChangePoints[rayI - 1].size();
             }
         }

         // Write points
         rayVtkFile
             << "DATASET POLYDATA" << nl
             << "POINTS " << nRayPoints << " double" << endl;

         // Add ray points
         forAll(beamDirectionChangePoints, rayI)
         {
             forAll(beamDirectionChangePoints[rayI], i)
             {
                 rayVtkFile
                     << beamDirectionChangePoints[rayI][i].x() << " "
                     << beamDirectionChangePoints[rayI][i].y() << " "
                     << beamDirectionChangePoints[rayI][i].z() << endl;
             }
         }

         // Count the number of lines
         label nRayLines = 0;
         forAll(beamDirectionChangePoints, rayI)
         {
             // Note: we must add 1 as the VTK format requires it
             nRayLines += beamDirectionChangePoints[rayI].size() + 1;
         }

         // Write lines
         rayVtkFile
             << "LINES " << beamDirectionChangePoints.size() << " " << nRayLines
             << endl;

         forAll(beamDirectionChangePoints, rayI)
         {
             // Write the number of points in the line
             rayVtkFile
                 << beamDirectionChangePoints[rayI].size();

             // Write indices of points
             forAll(beamDirectionChangePoints[rayI], i)
             {
                 rayVtkFile
                     << " " << pointIdOffset[rayI] + i;
             }

             rayVtkFile
                 << endl;
         }
     }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
