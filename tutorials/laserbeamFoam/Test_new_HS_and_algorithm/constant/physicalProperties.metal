/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  10
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      physicalProperties.water;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

viscosityModel  constant;

nu              5e-07;

rho             8000;

    poly_kappa   (10 0.015 0 0 0 0 0 0);
    //poly_cp   (244.8 9.587e-1 -3.77e-4 6.5e-8 -4.14e-12 0 0 0);
    poly_cp   (520 0.075 0 0 0 0 0 0);
    
    elec_resistivity	1e-6;
    
	Tsolidus 1658;
	Tliquidus 1723;
    LatentHeat 2.7e5;
    beta    5.0e-6;


// ************************************************************************* //
