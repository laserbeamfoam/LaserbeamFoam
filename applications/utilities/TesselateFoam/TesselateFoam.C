
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.
    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.
    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.
Application
    VoronoiFoam
Description
    Utility to Set Phase Fields Using Voronoi Tesselation
Author
	Tom Flint
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// #include <random>
// #include <vector>
// #include <algorithm>
#include "./Voronoi/src/voro++.hh"

double rnd() {return double(rand())/RAND_MAX;}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{


#   include "setRootCase.H"
#   include "createTime.H"



#   include "createMesh.H"


    volScalarField xi
    (
        IOobject
        (
            "xi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
	    dimensionedScalar("xi",dimensionSet(0, 0, 0, 0, 0),0.0),
	    zeroGradientFvPatchScalarField::typeName
    );


    IOdictionary PhaseFieldProperties
    (
        IOobject
        (
            "PhaseFieldProperties",    // dictionary name
            runTime.constant(),     // dict is found in "constant"
            mesh,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE      // dict is only read by the solver
        )
    );


//label N_Ori(readLabel(PhaseFieldProperties.lookup("N_Phases")));
label N_Seeds(readLabel(PhaseFieldProperties.lookup("N_Seeds")));

label N_Ori = N_Seeds+1;

Info<<"Number of SEEDS: "<<N_Seeds<<endl;



PtrList<volScalarField> PopBal(N_Ori);

    for(label count = 0; count < N_Ori; count++)
{

    word name_ni ("n." + name(count));

    PopBal.set
    (
        count,
        new volScalarField 
    (
            IOobject
            (
                name_ni,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
	    dimensionedScalar("ni",dimensionSet(0, 0, 0, 0, 0),0.0),
	    zeroGradientFvPatchScalarField::typeName
        )
    );
}

scalar Xmin(readScalar(PhaseFieldProperties.lookup("Xmin")));
scalar Xmax(readScalar(PhaseFieldProperties.lookup("Xmax")));
scalar Ymin(readScalar(PhaseFieldProperties.lookup("Ymin")));
scalar Ymax(readScalar(PhaseFieldProperties.lookup("Ymax")));
scalar Zmin(readScalar(PhaseFieldProperties.lookup("Zmin")));
scalar Zmax(readScalar(PhaseFieldProperties.lookup("Zmax")));


double X0 = Xmin;//-0.05e-3;
double X1 = Xmax;
double Y0 = Ymin;
double Y1 = Ymax;//2.0e-3;
double Z0 = Zmin;
double Z1 = Zmax;

int particles=N_Seeds;


            voro::container con(X0,X1,Y0,Y1,Z0,Z1,10,10,1,
            false,false,false,8);

            
int K;
    // Randomly add particles into the container
    for(int i=0;i<particles;i++) {
        double x=X0+rnd()*(X1-X0);
        double y=Y0+rnd()*(Y1-Y0);
        double z=Z0+rnd()*(Z1-Z0);
        con.put(i,x,y,z);
    }


double r,rx,ry,rz;

forAll (mesh.C(), celli) {

if(con.find_voronoi_cell(mesh.C()[celli].component(0),mesh.C()[celli].component(1),mesh.C()[celli].component(2),rx,ry,rz,K)){
// GRID[x[0]][x[1]][x[2]][K/N_coincident_phases]=1.0;
if(PopBal[N_Seeds][celli]<1.0-SMALL && xi[celli]>1.0-SMALL){
PopBal[K][celli]=1.0;
}
}

}



for(label count = 0; count < N_Ori; count++)//find the number of non zero phases at each point in the domain
{
PopBal[count].write();
}
/*
double r,rx,ry,rz;
    for (int i=0; i<nodes(GRID); i++) {
vector<int> x = position(GRID,i);
GRID_MELTTRACK(i)[0]=0.0;
if(con.find_voronoi_cell(x[0],x[1],x[2],rx,ry,rz,K)){
GRID[x[0]][x[1]][x[2]][K/N_coincident_phases]=1.0;
}
}
*/

// Populate grainNum field
volScalarField grainNum
(
    IOobject
    (
        "grainNum",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("grainNum",dimensionSet(0, 0, 0, 0, 0),-1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField maxNiVal
(
    IOobject
    (
        "maxNiVal",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("maxNiVal",dimensionSet(0, 0, 0, 0, 0),0.0),
    zeroGradientFvPatchScalarField::typeName
);
scalar grainNumThreshold = 0.5;
forAll(grainNum,i){
    forAll(PopBal,j){
        if(
            (PopBal[j][i] > grainNumThreshold)
        &&  (PopBal[j][i] > maxNiVal[i])
        )
        {
            maxNiVal[i] = i;
            grainNum[i] = j;
        }   
    }
}
grainNum.write();


    return 0;
}


#include "./Voronoi/src/voro++.cc"

// ************************************************************************* //