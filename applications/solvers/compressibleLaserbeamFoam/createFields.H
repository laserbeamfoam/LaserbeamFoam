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

// EOS density sum_i alpha_i*rho_i(p,T): pEqn compressibility, sound speed, viscosity only
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

// floor for the ledger density in empty cells (same value as rhoMinEOS)
dimensionedScalar rhoMMin
(
    "rhoMMin",
    dimDensity,
    mixture.getOrDefault<scalar>("rhoMMin", 1e-3)
);

// inertia and heat capacity never below the EOS values the pressure equation assumes (ledgerFloorEOS false = pure ledger)
const bool ledgerFloorEOS(mixture.getOrDefault<bool>("ledgerFloorEOS", false));
Info<< "ledgerFloorEOS = "
    << (
           ledgerFloorEOS
         ? "true: rhoM, rhoCv = max(ledger, EOS), UEqn in rhoM*DU/Dt form"
         : "false: pure ledger"
       )
    << endl;

// transported density sum_i c_i, the mass rhoPhi moves (updateRhoM.H); rebuilt from c_i on restart
volScalarField rhoM
(
    IOobject
    (
        "rhoM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    max(mixture.rhoM(), rhoMMin)
);
if (ledgerFloorEOS)
{
    rhoM = max(rhoM, rho);
}
rhoM.oldTime();     // oldTime storage for fvm::ddt(rhoM, U)

mesh.setFluxRequired(p_rgh.name());


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


// Construct compressible turbulence model
// EOS rho on purpose: laminar rho*nuEff == mu only with the density nu was built from
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


       volScalarField vDot // phase-change dilatation the ledger applied [1/s]; its cells get full ddtCorr memory (pEqn.H)
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



        // latent heat of vaporisation for TEqn.H [W/m^3]: Lvap*r at the anchors plus the p*dV and sensible
        // heat of the vapour moved (r > 0 evaporates). Rebuilt each step by mixture.solve(), so NO_READ
        // (a restart would otherwise read an old-dimension mass_dot and fail the dimension check)
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


        // volume closure chi*(s - 1)/dt [1/s], set after TEqn by volumeClosure.H; the pEqn's volume source
        volScalarField closureSrc
    (
        IOobject
        (
            "closureSrc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("closureSrc", dimensionSet(0,0,-1,0,0), 0.0)
    );


        // phase-change dilatation before the per-step volume cap [1/s], dt control only
        volScalarField vDotDt
    (
        IOobject
        (
            "vDotDt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("vDotDt", dimensionSet(0, 0, -1, 0, 0), 0.0)
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

// latent heat of fusion per volume, rhoLatentHeat = sum_{i condensed} c_i*L_i [J/m^3]: the volume-additive
// form of laserbeamFoam's rho*L_mass (LatentHeat*rho under-counted an interface cell by about the factor
// condensate). Filled in update.H once per step
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

// rhoCv = sum_i c_i*Cv_i [J/m^3/K]: the heat capacity per unit
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

// CvMass = rhoCv/rhoM = sum_i Y_i*Cv_i [J/kg/K], the mass-weighted mixture Cv (Y_i = c_i/rhoM);
// written so it can be checked: ~Cv_metal in metal, ~Cv_gas in gas
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


// mushy-zone constant C [kg/m^3/s] and guard b of the Darcy drag (UEqn.H), optional in fvSolution
// MELTING. Carman-Kozeny: C = 180*mu/lambda2^2 for dendrite arm spacing lambda2 (REVIEW.md)
const dimensionedScalar DarcyConstantlarge
(
    "DarcyConstantlarge", dimensionSet(1,-3,-1,0,0),
    mesh.solutionDict().subDict("MELTING").getOrDefault<scalar>("DarcyConstant", 1.0e6)
);

const dimensionedScalar DarcyConstantsmall
(
    "DarcyConstantsmall", dimless,
    mesh.solutionDict().subDict("MELTING").getOrDefault<scalar>("DarcyConstantSmall", 1.0e-12)
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
#include "createMRF.H"
#include "createK.H"

// ---- Stefan flow of each step's condensation (stefanFlow.H, UEqn.H, pEqn.H) ----------------------
// psiStefan: its potential; fixed (0) where p_rgh is fixed, so the flow leaves through the open
// boundaries, zero gradient elsewhere (constraint patches keep their type)
wordList psiStefanTypes(p_rgh.boundaryField().types());
forAll(psiStefanTypes, patchi)
{
    if (p_rgh.boundaryField()[patchi].fixesValue())
    {
        psiStefanTypes[patchi] = fixedValueFvPatchScalarField::typeName;
    }
    else
    {
        psiStefanTypes[patchi] = zeroGradientFvPatchScalarField::typeName;
    }
}

volScalarField psiStefan
(
    IOobject
    (
        "psiStefan",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimArea/dimTime, Zero),
    psiStefanTypes
);
mesh.setFluxRequired(psiStefan.name());

// its cell velocity (a donor carries its gas side's), weighted by rho_gas/rho for the viscous term
volVectorField Upc
(
    IOobject
    (
        "Upc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector(dimVelocity, Zero),
    zeroGradientFvPatchVectorField::typeName
);

// the Stefan flux the linear interpolation of Upc misses, on the donors' faces
surfaceScalarField phiJump
(
    IOobject
    (
        "phiJump",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar(dimVolume/dimTime, Zero)
);
phiJump.setOriented();

// the same per unit face area, kept one step for the ddtCorr memory (pEqn.H); AMR maps it face by
// face, so a split face keeps the velocity, not the parent's whole flux. Written, so a restart
// continues with the jump the old flux carries
surfaceScalarField uJump
(
    IOobject
    (
        "uJump",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimVelocity, Zero)
);
uJump.setOriented();               // a normal velocity: its sign follows the face (AMR flips)
uJump.oldTime().setOriented();
