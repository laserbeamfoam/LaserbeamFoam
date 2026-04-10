dimensionedScalar deltaN=1e-8/pow(average(mesh.V()), 1.0/3.0);


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Constructing multiphaseMixtureThermo\n" << endl;
multiphaseMixtureThermo mixture(U, phi);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mixture.rho()
);

dimensionedScalar pMin("pMin", dimPressure, mixture);

mesh.setFluxRequired(p_rgh.name());


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


// Construct compressible turbulence model
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        mixture.rhoPhi(),
        mixture
    )
);


       volScalarField vDot //kg per m^3 per s
    (
        IOobject
        (
            "vDot",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
       dimensionedScalar("vDot",dimensionSet(0,0,-1,0,0),0.0)
    );

    

        volScalarField mass_dot //
    (
        IOobject
        (
            "mass_dot",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("mass_dot", dimensionSet(1, -3, -1, 1, 0), 0.0)
    //    dimensionedScalar("mass_dot",dimensionSet(1,-1,-3,0,0),0.0)
    );


     volScalarField alphasum
    (
        IOobject
        (
            "alphasum",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
       dimensionedScalar("alphasum",dimensionSet(0,0,0,0,0),0.0)
    );

        volScalarField condensate
    (
        IOobject
        (
            "condensate",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("condensate",dimless,0.0)
    );

            volScalarField condensateFiltered
    (
        IOobject
        (
            "condensateFiltered",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("condensateFiltered",dimless,0.0)
    );
    condensateFiltered = max(min(condensate, scalar(1)), scalar(0));
    condensateFiltered.correctBoundaryConditions();

    //         volScalarField condensatesmooth
    // (
    //     IOobject
    //     (
    //         "condensatesmooth",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //     dimensionedScalar("condensatesmooth",dimless,0.0)
    // );
    // condensatesmooth=fvc::average(condensate);


    volScalarField metalaverage
(
    IOobject
    (
        "metalaverage",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("metalaverage",dimensionSet(0,0,0,0,0),0.0)
);
metalaverage = fvc::average(condensate);


    // volScalarField refine_test //kg per m^3 per s
    // (
    //     IOobject
    //     (
    //         "refine_test",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::READ_IF_PRESENT,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //    dimensionedScalar("refine_test",dimensionSet(0,0,0,0,0),0.0)
    // );

volVectorField n_filtered
(
    IOobject
    (
        "n_filtered",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("n_filtered",dimensionSet(0, 0, 0, 0, 0),vector::zero)
);









volScalarField TSolidus
    (
        IOobject
        (
            "TSolidus",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("TSolidus",dimensionSet(0, 0, 0, 1, 0),1600)
    );

          volScalarField TLiquidus
    (
        IOobject
        (
            "TLiquidus",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("TLiquidus",dimensionSet(0, 0, 0, 1, 0),1800.0)
    );

    volScalarField LatentHeat
    (
        IOobject
        (
            "LatentHeat",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("LatentHeat",dimensionSet(0, 2, -2, 0, 0),2.7e5)
        //(cp*(Temperature-((fluid.Tliquidus()+fluid.Tsolidus())/(scalar(2)))) + delta1*fluid.Latentheat())
    );


        volScalarField electrical_resistivity//electrical resistivity of condensate only
    (
        IOobject
        (
            "electrical_resistivity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("electrical_resistivity",dimensionSet(1,3,-3,0,0),0.0)
    );







autoPtr<mthdModel> mthd(nullptr);

if (mthdModel::electromagneticsEnabled(mesh, "thermophysicalProperties"))
{
    mthd.reset(new mthdModel(mesh, "thermophysicalProperties"));
}


const dimensionedScalar DarcyConstantlarge
(
    "DarcyConstantlarge", dimensionSet(1,-3,-1,0,0), 1.0e6
);

const dimensionedScalar DarcyConstantsmall
(
    "DarcyConstantsmall", dimless, 1.0e-12
);






// Create laser heat source
laserHeatSource laser(mesh);


#include "update.H"



// TSolidus*=0.0;
// TLiquidus*=0.0;
// LatentHeat*=0.0;

// forAllConstIter(PtrDictionary<phaseModel>, mixture.phases(), iter)
//     {
//         const phaseModel& phase1 = iter();
//         const volScalarField& alpha1 = phase1;

//         // mesh.setFluxRequired(phase1.name());

//         IOdictionary phasedictionary
//         (
//             IOobject
//             (
//             "thermophysicalProperties." + phase1.name(),
//             mesh.time().constant(),
//             mesh,
//             IOobject::MUST_READ_IF_MODIFIED
//             )
//         );
        
//         // dimensionedScalar molweight("molweight",dimless,specieDict.lookup("molWeight"));
//         // phasedictionary.lookup("TSolidus");
//         // dimensionedScalar phaseTSolidus("phaseTSolidus",dimensionSet(0, 0, 0, 1, 0),phasedictionary.lookup("TSolidus"));
//         dimensionedScalar phaseTSolidus("phaseTSolidus",dimensionSet(0, 0, 0, 1, 0),readScalar(phasedictionary.lookup("TSolidus")));
//         dimensionedScalar phaseTLiquidus("phaseTLiquidus",dimensionSet(0, 0, 0, 1, 0),readScalar(phasedictionary.lookup("TLiquidus")));
//         dimensionedScalar phaseLH("phaseLH",dimensionSet(0, 2, -2, 0, 0),readScalar(phasedictionary.lookup("LatentHeat")));

//         TSolidus+=alpha1*phaseTSolidus;
//         TLiquidus+=alpha1*phaseTLiquidus;
//         LatentHeat+=alpha1*phaseLH;

//     }





//     TSolidus=mixture.Tsolidus();
// TLiquidus=mixture.Tliquidus();

        volScalarField DC
    (
        IOobject
        (
            "DC",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("DC",dimensionSet(1,-3,-1,0,0),1.0e14)
    );

        volScalarField epsilon1
    (
        IOobject
        (
            "epsilon1",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
       dimensionedScalar("epsilon1",dimensionSet(0,0,0,0,0),0.0)
    );
    
    epsilon1=max(min((mixture.T()-TSolidus)/(TLiquidus-TSolidus),scalar(1)),scalar(0));
    epsilon1.correctBoundaryConditions();


         volScalarField epsilon1mask
    (
        IOobject
        (
            "epsilon1mask",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
       dimensionedScalar("epsilon1mask",dimensionSet(0,0,0,0,0),0.0)
    );



volScalarField TRHS
(
    IOobject
    (
        "TRHS",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedScalar("TRHS",dimensionSet(1, -1, -3, 0, 0),scalar(0.0))
);

volScalarField Tcorr
    (
        IOobject
        (
            "Tcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("Tcorr",dimensionSet(0, 0, 0, 1, 0),0.0)
    );

                volScalarField DilationError
    (
        IOobject
        (
            "DilationError",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("DilationError",dimensionSet(0,0,-1,0,0),0.0)
    );

    

volVectorField Marangoni
(
    IOobject
    (
        "Marangoni",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("Marangoni",dimensionSet(1, -2, -2, 0, 0),vector::zero)
);

const volVectorField zeroJcrossB
(
    IOobject
    (
        "zeroJcrossB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zeroJcrossB", dimensionSet(1, -2, -2, 0, 0), vector::zero)
);

const volScalarField zeroJoulesource
(
    IOobject
    (
        "zeroJoulesource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroJoulesource", dimensionSet(1, -1, -3, 0, 0), 0.0)
);


#include "createMRF.H"
#include "createK.H"
