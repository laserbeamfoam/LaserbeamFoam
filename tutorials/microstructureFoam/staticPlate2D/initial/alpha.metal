FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      alpha.metal;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    atmosphere
    {
        type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;
    }

    frontAndBack
    {
        type            empty;
    }

    rightWall
    {
        type            zeroGradient;
    }

    leftWall
    {
        type            zeroGradient;
    }

    lowerWall
    {
        type            zeroGradient;
    }
}


// ************************************************************************* //
