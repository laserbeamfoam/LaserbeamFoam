/*--------------------------------*- C++ -*----------------------------------*\
  =========
  \\      /  F ield           OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration       Website:  www.openfoam.com
    \\  /    A nd             Version:  v2506
     \\/     M anipulation
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       volScalarField;
    location    "0";
    object      alpha.gas;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    inlet
    {
        type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;
    }

    outlet
    {
        type            zeroGradient;
    }

    lowerWall
    {
        type            zeroGradient;
    }

    upperWall
    {
        type            zeroGradient;
    }

    defaultFaces
    {
        type            empty;
    }
}


// ************************************************************************* //
