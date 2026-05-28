/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
\*---------------------------------------------------------------------------*/

#include "movingConeInjection.H"
#include "basicKinematicCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef basicKinematicCloud::kinematicCloudType kinematicCloudType;
    
    defineTemplateTypeNameAndDebug(movingConeInjection<kinematicCloudType>, 0);
    
    InjectionModel<kinematicCloudType>::
        adddictionaryConstructorToTable<movingConeInjection<kinematicCloudType>>
        addmovingConeInjectionkinematicCloudTypeConstructorToTable_;
}

// ************************************************************************* //