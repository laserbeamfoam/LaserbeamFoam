/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  3.0.1                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      alpha;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{

	// #includeEtc "caseDicts/setConstraintTypes"
    crack
    {
        type            zeroGradient;
    }
	
	// 2. The interface created by zoning (Side A)
    metalZone_to_airZone
    {
        type            zeroGradient; 
    }
	
    // 3. The interface created by zoning (Side B)
    airZone_to_metalZone
    {
        type            zeroGradient;
    }
    
	frontAndBack{
        type            zeroGradient;
	}
	
    lowerWall
    {
        type            zeroGradient;
    }
	 
    atmosphere
    {
        type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;
    }
    rightWall
    {
        type            zeroGradient;
    }
    leftWall
    {
        type            zeroGradient;
    }
	
    defaultFaces
    {
        type            empty;
    }
	 
	".*"
    {
        type            processor;
        value           uniform 0;
    }
}
// ************************************************************************* //
