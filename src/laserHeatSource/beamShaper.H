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
    outward-facing unit normal of (0 1 ).

    Laser rays are written as VTKs for visualisation.

Author
    Thomas Flint, The University of Manchester.
    Philip Cardiff, University College Dublin.

SourceFiles
    laserHeatSource.C

\*---------------------------------------------------------------------------*/

#ifndef beamShaper_H
#define beamShaper_H

#include "volFields.H"
#include "dictionary.H"
#include "scalar.H"
#include "List.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{ 
/*---------------------------------------------------------------------------*\
                        Class laserHeatSource Declaration
\*---------------------------------------------------------------------------*/

class BeamShape
{
    /***********************************\
     * @brief: Class storing/producing * 
     * the matrix address list data:   *

                 ↑ 
                +z 
            ← -x o +x →    
                -z 
                 ↓ 

       The address book is 33x33:

        ↑  ↑  ↑  ↑  ↑  ↑  ↑  ↑  ↑  
       -,+ ·  ·  · 0,+ ·  ·  · +,+     
        ·  ·  ·  ·  ·  ·  ·  ·  ·      
        ·  ·  ·  ·  ·  ·  ·  ·  ·      
        ·  ·  ·  ·  ·  ·  ·  ·  ·      
       -,0 ·  ·  ·  ●  ·  ·  · +,0     
        ·  ·  ·  ·  ·  ·  ·  ·  ·      
        ·  ·  ·  ·  ·  ·  ·  ·  ·      
        ·  ·  ·  ·  ·  ·  ·  ·  ·       
       -,- ·  ·  · 0,- ·  ·  · +,- 
     
    \************************************/

private:

    /*********************************\
     * Beam shape file and addresses *
    \*********************************/
    bool debug = false;                           // Debug flag
    fileName shapeFile_;                         // Filename for the beam shape
    fileName scheduleFile_;                      // Per-profile beam-shape schedule file
    const fvMesh& mesh_; // If you also need mesh to be accessible later on
    const dictionary& laserProperties_; // Add this line to store LaserProperties
    DynamicList<word> beamShapesList_;                  // This is read in from the laser properties file
    DynamicList<List<Pair<label>>> beamShapeLibrary_;   // List of beam shape addresses, organized by shape
    DynamicList<Pair<label>> beamShapeAddresses_;       // List of beam shape addresses, organized by time, ++switch state
    /*********************************\
     * Beam shape timing information *
    \*********************************/
    Pair<label> dummyAddress_;       // Dummy address for the beam shape
    List<scalar> beamShapeTiming_;   // List of times for each beam shape
    List<label> beamShapeIndex_;     // List of indices for each beam shape
    scalar addressFrequency_;        // Frequency of address change
    List<scalar> addressDuration_;   // Duration of address change
    List<scalar> shapeFrequency_;    // beamShapeFile Frequency how many times the beam shape iterated Hz
    List<scalar> shapeDuration_;     // Duration of each beam shape exposure
    label shapeRepetitions_;        // Number of times the shape is repeated in shapeDuration_
    scalar cycleDuration_;           // total time of all shapes in cycle
    scalar addressInterchange_       = 1E-8;  /* Time between address changes, 
                                                        should this be computed depending on:
                                                        - address distance? r2?
                                                        - OPA controller time? - how to measure this */ 
    
    /*********************************\
     * Storage array static sizes    *
    \*********************************/
    const label dim_ = 550;                         // Side length of beam shape in pixels
    const label dim_grd = 33;                       // Side length of address grid in pixels
    const label dim_inp = 203;                      // Side length of pirmitive shape in pixels
    const label dim_exp = 1429;                     // Side length of expanded beam shape in pixels
    dimensionedScalar pixelSize_;                   // Pixel size, same in x and z in metres
    dimensionedScalar pixelArea_;                   // Pixel area in m^2
    
    List<scalar> beamBox_;                          // storage for the generator distribution
    List<scalar> addressShape_;                     // storage for the current address shape 
    List<scalar> beamShape_;                        // storage for the current beam shape
    List<scalar> GaussianField_;                    // storage for the Gaussian filter
    const scalar sigma = 550/6.6;                   /* Sigma value for Gaussian filter (550/6.6), 
                                                        should we use a multivariate Gaussian? */
    HashSet<fileName> writtenBeamShapeFiles_;        // Track CSV files already written this run
    /***********************************\
     * Member funcs: set,get,mod,Gauss *
    \***********************************/
    
    scalar getValue(List<scalar>& arrayData, const Pair<label>& coords, label dim);
    scalar getValue(List<scalar>& arrayData, label index);
    void setValue(List<scalar>& arrayData, label dim, const Pair<label>& coords, const scalar value);
    void setValue(List<scalar>& arrayData, label index, scalar value);

    void generateGaussianFilter(scalar sigma);
    void applyGaussianFilter(List<scalar>& arrayData, bool inverse);
    void generateExpandedBeamShape(const label num_copies);
    void normalizeField(List<scalar>& arrayData);
    void generateMeanImage();
    void printBeamShape(const fileName& filename, List<scalar>& arrayData, label dim);
    
    List<word> readShapeNames(const dictionary& dict);
    List<Pair<label>> readAddressesByName(const dictionary& dict, const word& listName);
    
    inline label index(const Pair<label>& coords, const label dim);
    Pair<label> coords(const label index, const label dim);

    void beamShapeFromAddress(const Pair<label>& beamShapeAddress);
    label currentAddressIndex(scalar currentTime);
    Pair<label> currentBeamShapeAddress(scalar currentTime);
    void shapeOverInterval(scalar currentTime, scalar timeStep);
    void beamShapeScheduler(const dictionary& beamShapeDict);
    void generateMeanImage(List<label> addressIndices, List<scalar> addressWeights);
    void beamShapeGenerator();

    // Add this declaration to allow updating internal fields to match the mesh
    // void updateMeshFields(const fvMesh& mesh);
    
public:
    BeamShape
    (
        const Foam::fileName& shapeFile,
        const Foam::fileName& scheduleFile,
        const fvMesh& mesh,
        const dictionary& laserProperties
    );
    ~BeamShape();

    dimensionedScalar imageLength_;                 // Length of the image in x and z (550*pixelSize) in metres
    scalar interpolateValue(scalar x_target, scalar z_target, scalar localCellSize, label dim); 
    void getBeamShape(scalar currentTime, scalar timeStep);
    // void updateMeshFields(const fvMesh& mesh);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
