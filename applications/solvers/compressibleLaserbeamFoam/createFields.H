// Check for restarting solver
// bRestartFirstLoop is used within `main()` also!
bool bRestartFirstLoop = (runTime.timeIndex() > 0);
Info<< "Runtime timeIndex = " << runTime.timeIndex() << endl;
if (bRestartFirstLoop)
{
    // solver was restarted!
    Info<< "\nSolver restart detected. Preventing correctPhi\'s on the first "
        << "pass/loop, and restoring `.oldTime()` values as required by "
        << "`fvc::ddt` methods (does not support fvc backwards scheme!).\n"
        << endl;
}

// Local debugging ------------------------------------------------------------
// Add debug entry for the solver (mainly to enable verbose fields)
Switch debugFlag
(
    IOdictionary
    (
        IOobject
        (
            "controlDict",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ).lookupOrDefault<Switch>("debugCompressibleLaserbeamFoam", false)
);

// Disabled temporarily:
// autoPtr<volScalarField> rhoCopyPreUpdatePtr; // reset after rho is constructed

// ----------------------------------------------------------------------------

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
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
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



        // Latent heat of vaporisation/condensation sink for TEqn.H:
        //   mass_dot = sum_pairs Lvap*r   [W/m^3]   (r > 0 evaporates)
        // Rebuilt from scratch by mixture.solve() every time step, so it is
        // NO_READ: the energy equation is now in rho*Cv form and this term
        // enters it directly in W/m^3 (it used to be rCv*Lvap*r in
        // [K kg/m^3/s]); a restart that read an old-dimension mass_dot file
        // would otherwise die on a dimension check in solveAlphas.
        volScalarField mass_dot
    (
        IOobject
        (
            "mass_dot",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("mass_dot", dimensionSet(1, -1, -3, 0, 0), 0.0)
    );


        volScalarField vDotP // -d(vDot)/dp >= 0: implicit phase-change pressure-coupling coefficient for pEqn,  [1/(s.Pa)] = dimensionSet(-1,1,1,0,0).  Filled by mixture.solve() once per time step;
                            
    (
        IOobject
        (
            "vDotP",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("vDotP", dimensionSet(-1, 1, 1, 0, 0), 0.0)
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
            IOobject::READ_IF_PRESENT,
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




volScalarField refineIndicator
(
    IOobject
    (
        "refineIndicator",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("refineIndicator", dimless, 0.0)
);




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
            IOobject::READ_IF_PRESENT,
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
            IOobject::READ_IF_PRESENT,
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
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("LatentHeat",dimensionSet(0, 2, -2, 0, 0),2.7e5)
        //(cp*(Temperature-((fluid.Tliquidus()+fluid.Tsolidus())/(scalar(2)))) + delta1*fluid.Latentheat())
    );


// ---- mass-weighted energy mixing (see updateRhoCv.H / update.H / TEqn.H) --

// Latent heat of fusion per unit VOLUME of the cell:
//   rhoLatentHeat = sum_{i condensed} alpha_i*rho_i*L_i   [J/m^3]
// This is the volume-additive quantity (the compressible analogue of
// laserbeamFoam's massFraction-weighted LatentHeat: rho*L_mass). The old
// product LatentHeat*rho = (sum_cond alpha_i L_i)*(sum_all alpha_j rho_j)
// under-counts the fusion latent heat of an interface cell by ~a factor
// condensate (10x in a 10% metal cell). Filled in update.H once per step
// (alphas are fixed for the step; the condensed rho_i barely move).
volScalarField rhoLatentHeat
(
    IOobject
    (
        "rhoLatentHeat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoLatentHeat", dimEnergy/dimVolume, 0.0)
);

// rhoCv = sum_i alpha_i*rho_i*Cv_i [J/m^3/K]: the heat capacity per unit
// volume that multiplies DT/Dt in TEqn.H (rho*Cv form, like laserbeamFoam's
// rhoCp). Needs an oldTime for fvm::ddt(rhoCv, T) / fvc::ddt(rhoCv).
volScalarField rhoCv
(
    IOobject
    (
        "rhoCv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoCv", dimEnergy/dimVolume/dimTemperature, 0.0)
);

// condensed-only part, for the source-based melt-fraction correction
// coefficient epsilonRel*rhoCvCond/rhoLatentHeat (= epsilonRel*Cv/L of the
// condensed material, mass weighted)
volScalarField rhoCvCond
(
    IOobject
    (
        "rhoCvCond",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoCvCond", dimEnergy/dimVolume/dimTemperature, 0.0)
);

// CvMass = rhoCv/rho = sum_i Y_i Cv_i [J/kg/K], the mass-weighted mixture
// Cv. Written so you can eyeball it: ~Cv_metal in metal, ~Cv_gas in gas.
// Its face interpolate times rhoPhi is the Cv-weighted mass flux of TEqn.H.
volScalarField CvMass
(
    IOobject
    (
        "CvMass",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("CvMass", dimEnergy/dimMass/dimTemperature, 0.0)
);


        volScalarField electrical_resistivity//electrical resistivity of condensate only
    (
        IOobject
        (
            "electrical_resistivity",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
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

// seed rhoCv/CvMass and enable their oldTime storage (first assignment in
// a new time step then shifts the end-of-previous-step values into
// .oldTime(), exactly as laserbeamFoam does with rhoCp.oldTime())
#include "updateRhoCv.H"
rhoCv.oldTime();
CvMass.oldTime();



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

    if (!epsilon1.headerOk())
    {
        Info<< "\nepsilon1: assuming fresh start, calculating the field...\n"
            << endl;

        // TSolidus/TLiquidus are condensed-only 
        epsilon1 = max
        (
            min
            (
                (mixture.T() - TSolidus)
               /max(TLiquidus - TSolidus, dimensionedScalar("dTfloor", dimTemperature, SMALL)),
                scalar(1)
            ),
            scalar(0)
        );
    }


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
        IOobject::NO_READ, //READ_IF_PRESENT,
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
            IOobject::READ_IF_PRESENT,
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
            IOobject::READ_IF_PRESENT,
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

// --- Perform oldTime restoration --------------------------------------------
if (bRestartFirstLoop)
{
    Info<< "Restoring oldTime values." << endl;
    int errorCount = 0;

    mixture.restoreOldTimeValues();

    volVectorField UOldTime
    (
        IOobject
        (
            IOobject::groupName("U", "oldTime"),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        U
    );

    if (UOldTime.headerOk())
    {
        U.oldTime() = UOldTime;
    }
    else
    {
        Info<<"Restart detected but U.oldTime not found!" << endl;
        errorCount++;
    }


    surfaceScalarField phiOldTime
    (
        IOobject
        (
            IOobject::groupName("phi", "oldTime"),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(UOldTime)
    );

    if (phiOldTime.headerOk())
    {
        phi.oldTime() = phiOldTime;
    }
    else
    {
        Info<< "Restart detected but phi.oldTime not found; reconstructing "
            << "it from U.oldTime." << endl;
        phi.oldTime() = phiOldTime;
        errorCount++;
    }


    volScalarField rhoOldTime
    (
        IOobject
        (
            IOobject::groupName("rho", "oldTime"),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        rho
    );

    if (rhoOldTime.headerOk())
    {
        rho.oldTime() = rhoOldTime;
    }
    else
    {
        Info<<"Restart detected but rho.oldTime not found!" << endl;
        errorCount++;
    }


    volScalarField LatentHeatOldTime
    (
        IOobject
        (
            IOobject::groupName("LatentHeat", "oldTime"),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        LatentHeat
    );

    if (LatentHeatOldTime.headerOk())
    {
        LatentHeat.oldTime() = LatentHeatOldTime;
    }
    else
    {
        Info<<"Restart detected but LatentHeat.oldTime not found!" << endl;
        errorCount++;
    }


    volScalarField epsilon1OldTime
    (
        IOobject
        (
            IOobject::groupName("epsilon1", "oldTime"),
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        epsilon1
    );

    if (epsilon1OldTime.headerOk())
    {
        epsilon1.oldTime() = epsilon1OldTime;
        epsilon1.correctBoundaryConditions();
    }
    else
    {
        Info<< "Restart detected but epsilon1.oldTime not found!" << endl;
        errorCount++;
    }

    if (errorCount)
    {
        Info<< "WARNING: " << errorCount << " fields' oldTime values were not "
            << "restored! This will cause issues within any fvc::ddt terms!"
            << "RECOMMENDED: please inspect your latest time folder!"
            << endl;
    }
    else
    {
        Info<< "Main: All relevant fields' oldTime values were restored using "
            << "file reads. All good!"
            << endl;
    }
}
// ----------------------------------------------------------------------------
// Debug fields related

// Disabled temporarily:
// if (debugFlag)
// {
//     rhoCopyPreUpdatePtr.reset
//     (
//         new volScalarField
//         (
//             IOobject
//             (
//                 "rhoCopyPreUpdate",
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             (rho)
//         )
//     );
// }
// ----------------------------------------------------------------------------

#include "createMRF.H"
#include "createK.H"