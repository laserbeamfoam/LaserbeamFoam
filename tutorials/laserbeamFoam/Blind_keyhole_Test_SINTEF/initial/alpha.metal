/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2506                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       volScalarField;
    object      alpha.metal;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{
	//#includeEtc "caseDicts/setConstraintTypes"

    atmosphere
    {
		//type            zeroGradient;
		
		type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;
        
    }
    
    lowerWall
    {
		type            zeroGradient;
    }    

    leftWall
    {
		type            zeroGradient;
    }
    
    rightWall
    {
		type            zeroGradient;
    }

	frontAndBack
	{
        type            zeroGradient;
	}

     
}


// ************************************************************************* //
