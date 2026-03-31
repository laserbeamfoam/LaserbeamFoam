#include "createRDeltaT.H"

const dimensionedScalar deltaN = 1e-8/pow(average(mesh.V()), 1.0/3.0);

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

Info<< "Reading field Temperature\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const word interfaceTrackingScheme
(
    transportProperties.lookup("interfaceTrackingScheme")
);

if
(
    interfaceTrackingScheme != "MULES"
 && interfaceTrackingScheme != "isoAdvector"
)
{
    FatalErrorIn("laserbeamFoam.C")
        << "Unknown interface tracking scheme "
        << interfaceTrackingScheme
        << ". Available options are 'MULES' and 'isoAdvector'"
        << exit(FatalError);
}

dictionary& transportPropertiesMetal = transportProperties.subDict("metal");
dictionary& transportPropertiesGas = transportProperties.subDict("gas");


const dimensionedScalar Tsolidus1
(
    "Tsolidus",
    dimensionSet(0, 0, 0, 1, 0),
    transportPropertiesMetal//.lookup("Tsolidus")
);
// Info<<Tsolidus1<<endl;
const dimensionedScalar Tsolidus2
(
    "Tsolidus",
    dimensionSet(0, 0, 0, 1, 0),
    transportPropertiesGas//.lookup("Tsolidus")
);
// Info<<Tsolidus2<<endl;

const dimensionedScalar Tliquidus1
(
    "Tliquidus",
    dimensionSet(0, 0, 0, 1, 0),
    transportPropertiesMetal//.lookup("Tliquidus")
);

const dimensionedScalar Tliquidus2
(
    "Tliquidus",
    dimensionSet(0, 0, 0, 1, 0),
    transportPropertiesGas//.lookup("Tliquidus")
);

const dimensionedScalar LatentHeat1
(
    "LatentHeat",
    dimensionSet(0, 2, -2, 0, 0),
    transportPropertiesMetal//.lookup("LatentHeat")
);

const dimensionedScalar LatentHeat2
(
    "LatentHeat",
    dimensionSet(0, 2, -2, 0, 0),
    transportPropertiesGas//.lookup("LatentHeat")
);

const dimensionedScalar beta1
(
    "beta",
    dimensionSet(0, 0, 0, -1, 0),
    transportPropertiesMetal//.lookup("beta")
);

const dimensionedScalar beta2
(
    "beta",
    dimensionSet(0, 0, 0, -1, 0),
    transportPropertiesGas//.lookup("beta")
);

const dimensionedScalar Marangoni_Constant
(
    "Marangoni_Constant",
    dimensionSet(1, 0, -2, -1, 0),
    transportProperties//.lookup("dsigmadT")
);

const dimensionedScalar p0
(
    "p0",
    dimensionSet(1, -1, -2, 0, 0, 0, 0),
    transportProperties//.lookup("p0")
); // atmospheric pressure

const dimensionedScalar Tvap
(
    "Tvap",
    dimensionSet(0, 0, 0, 1, 0),
    transportProperties//.lookup("Tvap")
);

const dimensionedScalar Mm
(
    "Mm",
    dimensionSet(1, 0, 0, 0, -1),
    transportProperties//.lookup("Mm")
);

const dimensionedScalar LatentHeatVap
(
    "LatentHeatVap",
    dimensionSet(0, 2, -2, 0, 0),
    transportProperties//.lookup("LatentHeatVap")
);

const dimensionedScalar R
(
    "R",
    dimensionSet(1, 2, -2, -1, -1, 0, 0),
        8.314
);

volScalarField em
(
    IOobject
    (
        "em",                       
        runTime.timeName(),         
        mesh,                       
        IOobject::NO_READ,          
        IOobject::NO_WRITE          
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0) 
);

const dimensionedScalar SB
(
    "SB",
    dimensionSet(1, 0, -3, -4, 0, 0, 0),
    scalar(5.67e-8)
); // Stefan-Boltzmann Constant

[[maybe_unused]]const scalar emL(readScalar(transportProperties.lookup("emL")));
[[maybe_unused]]const scalar emS(readScalar(transportProperties.lookup("emS")));
const dimensionedScalar TRef ("TRef",transportProperties);

const Polynomial<8> polykappa_m(transportPropertiesMetal.lookup("poly_kappa"));
const Polynomial<8> polykappa_g(transportPropertiesGas.lookup("poly_kappa"));

const Polynomial<8> polycp_m(transportPropertiesMetal.lookup("poly_cp"));
const Polynomial<8> polycp_g(transportPropertiesGas.lookup("poly_cp"));

volScalarField cp
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//change to no write when working
    ),
    mesh,
    dimensionedScalar("cp", dimensionSet(0, 2, -2, -1, 0), 0.0)
);

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//change to no write when working
    ),
    mesh,
    dimensionedScalar("kappa",dimensionSet(1, 1, -3, -1, 0),0.0)
);

#include "updateKappaCp.H"


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
    alpha1*Tsolidus1 + (1.0 - alpha1)*Tsolidus2
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
    alpha1*Tliquidus1 + (1.0 - alpha1)*Tliquidus2
);

volScalarField DC
(
    IOobject
    (
        "DC",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DC", dimensionSet(1,-3,-1,0,0), 1.0e14)
);

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
    fvc::average(alpha1)
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
    max(min((T - TSolidus)/(TLiquidus - TSolidus), 1.0), 0.0)
);

volVectorField n_filtered
(
    IOobject
    (
        "n_filtered",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("n_filtered", dimless, vector::zero)//,
    // zeroGradientFvPatchVectorField::typeName
);

volScalarField epsilon1mask
(
    IOobject
    (
        "epsilon1mask",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilon1mask",dimless,0.0)
);

volScalarField epsilon1mask2
(
    IOobject
    (
        "epsilon1mask2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilon1mask2", dimless, 0.0)
);

volVectorField gradT
(
    IOobject
    (
        "gradT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("gradT", dimensionSet(0, -1, 0, 1, 0), vector::zero)
);

volScalarField LatentHeat
(
    IOobject
    (
        "LatentHeat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*LatentHeat1 + (1.0 - alpha1)*LatentHeat2
);

volScalarField TRHS
(
    IOobject
    (
        "TRHS",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("TRHS", dimensionSet(1, -1, -3, 0, 0), 0.0)
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
    dimensionedScalar("Tcorr", dimensionSet(0, 0, 0, 1, 0), 0.0)
);

volScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*beta1 + (1.0-alpha1)*beta2
);

volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhok", dimless, 0.0)
);

volScalarField Num_divU
(
    IOobject
    (
        "Num_divU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Num_divU", dimensionSet(0,0,-1,0,0), 0.0)
);

volVectorField Marangoni
(
    IOobject
    (
        "Marangoni",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("Marangoni", dimensionSet(1, -2, -2, 0, 0), vector::zero)
);

// Recoil Pressure
volScalarField pVap
(
    IOobject
    (
        "pVap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
);

// Evaporation Cooling
volScalarField Qv
(
    IOobject
    (
        "Qv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1,0,-3,0,0,0,0), 0.0)
);


volScalarField damper
(
    IOobject
    (
        "damper",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("damper", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField thermalDamper
(
    IOobject
    (
        "thermalDamper",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("thermalDamper", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField alpha_filtered
(
    IOobject
    (
        "alpha_filtered",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha_filtered", dimensionSet(0, 0, 0, -0, 0), 0.0)
);

volScalarField alpha_smoothed
(
    IOobject
    (
        "alpha_smoothed",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha_smoothed", dimensionSet(0, 0, 0, -0, 0), 0.0)
);

// volVectorField gradepsilon1(fvc::grad(epsilon1));
// volScalarField e1temp(fvc::average(epsilon1));
// volScalarField e2temp(fvc::average(epsilon1mask));


volScalarField e1temp
(
    IOobject
    (
        "e1temp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("e1temp", dimensionSet(0, 0, 0, -0, 0), 0.0)
);
e1temp = fvc::average(epsilon1);


dimensionedScalar elec_resistivity_metal
(
    "elec_resistivity",
    dimensionSet(1, 3, -3, 0, 0),
    transportPropertiesMetal//.lookup("elec_resistivity")
);


// Electrical resistivity of condensate only
volScalarField electrical_resistivity
(
    IOobject
    (
        "electrical_resistivity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    elec_resistivity_metal
);

volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoCp", dimensionSet(1, -1, -2, -1, 0), 0.0)
);
rhoCp = (rho*cp);
// volScalarField rhoCp(rho*cp);
rhoCp.oldTime();


// surfaceScalarField rhophicp(fvc::interpolate(cp)*rhoPhi);

surfaceScalarField rhophicp
(
    IOobject
    (
        "rhophicp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhophicp", dimensionSet(1, 2, -3, -1, 0), 0.0)
);
rhophicp = fvc::interpolate(cp)*rhoPhi;











const dimensionedScalar DarcyConstantlarge
(
    "DarcyConstantlarge", dimensionSet(1,-3,-1,0,0), 1.0e6
);

const dimensionedScalar DarcyConstantsmall
(
    "DarcyConstantsmall", dimless, 1.0e-12
);



// volVectorField gradAlpha(fvc::grad(alpha1));

volVectorField gradAlpha
(
    IOobject
    (
        "gradAlpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("gradAlpha", dimensionSet(0, -1, 0, -0, 0), vector::zero)
);
gradAlpha = fvc::grad(alpha1);

// volVectorField nHatM(gradAlpha/(mag(gradAlpha) + deltaN));


volVectorField nHatM
(
    IOobject
    (
        "nHatM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("nHatM", dimensionSet(0, 0, 0, -0, 0), vector::zero)
);
nHatM = gradAlpha/(mag(gradAlpha) + deltaN);


// This variable is >=0 if the cell has been liquid
volScalarField meltHistory
(
    IOobject
    (
        "meltHistory",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("meltHistory", dimensionSet(0, 0, 0, -0, 0), 0.0)
);

// This variable is 1 if (alpha.metal >= 0.5 and epsilon1 >=0.5), 0 otherwise
volScalarField condition
(
    IOobject
    (
        "condition",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("condition", dimensionSet(0, 0, 0, -0, 0), 0.0)
);





//microstructure stuff

dimensionedScalar pi("pi", dimensionSet(0, 0, 0, 0, 0), constant::mathematical::pi);



    volScalarField GBs
    (
        IOobject
        (
            "GBs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("GBs",dimensionSet(0, 0, 0, 0, 0),0.0)
    );
        volScalarField Sum_Phases
    (
        IOobject
        (
            "Sum_Phases",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Sum_Phases",dimensionSet(0, 0, 0, 0, 0),0.0)
    );




        volScalarField Lg
    (
        IOobject
        (
            "Lg",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Lg",dimensionSet(-1, 1, 1, 0, 0),0.0)
    );

        volScalarField Lp
    (
        IOobject
        (
            "Lp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Lp",dimensionSet(-1, 1, 1, 0, 0),0.0)
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

#   include "setNucleationSites.H"
// This will set ni_fieldsNeeded to the number of nucleation sites assigned.

// label N_Ori(readLabel(PhaseFieldProperties.lookup("N_Phases")));
    // Need to get number of tesselated fields.
label N_Seeds(readLabel(PhaseFieldProperties.lookup("N_Seeds")));
label N_Ori = ni_fieldsNeeded + N_Seeds;

Info<<"Number of Phase Field Order Parameters: "<< N_Ori << nl
<< "\tSeeds: " << N_Seeds << nl
<< "\tNuc Sites: " << ni_fieldsNeeded << nl
<<endl;

Info<<"Number of Phase Field Order Parameters: "<<N_Ori<<endl;

// abort();


PtrList<volScalarField> PopBal(N_Ori);
PtrList<volScalarField> Sum_parts(N_Ori);

for(label count = 0; count < N_Ori; count++)
{

    word name_Sum_parts ("Sum_parts" + name(count));

    Sum_parts.set
    (
        count,
        new volScalarField
        (
            IOobject
            (
                name_Sum_parts,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Sum_parts",dimensionSet(0, 0, 0, 0, 0),0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

List<bool> niActive (PopBal.size(),true);

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
                IOobject::READ_IF_PRESENT,//IOobject::MUST_READ,
                // IOobject::AUTO_WRITE
                IOobject::NO_WRITE
            ),
            mesh,
	    dimensionedScalar("ni",dimensionSet(0, 0, 0, 0, 0),0.0),
	    zeroGradientFvPatchScalarField::typeName
        )
    );
}
//Phase Field Stuff for Non-Conserved Parameters
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


PtrList<volScalarField> PopBal_av(N_Ori);

for(label count = 0; count < N_Ori; count++)
{

    word name_ni_av ("n_av" + name(count));

    PopBal_av.set
    (
        count,
        new volScalarField
        (
            IOobject
            (
                name_ni_av,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("niav",dimensionSet(0, 0, 0, 0, 0),0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

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

scalar grainNumThreshold( readScalar(PhaseFieldProperties.lookup("grainNumThreshold")));

Switch write_ni_all( PhaseFieldProperties.lookupOrDefault<Switch>("write_ni_all",true) );
Switch write_ni_active( PhaseFieldProperties.lookupOrDefault<Switch>("write_ni_active",true) );


// Set grain num field to indicate which field occupies the cell
grainNum = dimensionedScalar(grainNum.dimensions(), -1); // Assume that no grains occupy any cell
maxNiVal = dimensionedScalar(maxNiVal.dimensions(), 0);

forAll(grainNum,i){
    forAll(PopBal,j){
        if(
            (PopBal[j][i] > grainNumThreshold)
        &&  (PopBal[j][i] > maxNiVal[i])
        )
        {
            maxNiVal[i] = PopBal[j][i];
            grainNum[i] = j;
        }   
    }
}


// Get max field number that currently exists
label maxNiIndex(-1);

forAll(PopBal,i){
    volScalarField& ni = PopBal[i];
    scalar maxVal = gMax(ni);
    Info << "PopBal[" <<i << "] max = " << maxVal << endl;
    if (maxVal > SMALL){ // field is non-zero
        // ni_active[i] = true; // turn on writing of field
        if (i > maxNiIndex)
        {
            maxNiIndex = i;
        }
    }
}

// the maxNiIndex that was found was specific to the processor
// Find the global maxNiIndex by maxOp
reduce(maxNiIndex,maxOp<label>());


volScalarField phi_PF
(
    IOobject
    (
        "phi_PF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("phi_PF",dimensionSet(0, 0, 0, 0, 0),0.0)
);



epsilon1=1.0-xi;

Random randObj2(clock::getTime());


double phi1[N_Ori]={0.0};

double phi2[N_Ori]={0.0};

double phi3[N_Ori]={0.0};

for(label count_j = 0; count_j < N_Ori; count_j++){

phi1[count_j]=randObj2.sample01<scalar>()*2.0*pi.value();
// phi2[count_j]=randObj2.scalar01()*2.0*pi.value();
// phi3[count_j]=randObj2.scalar01()*2.0*pi.value();

}




PtrList<volVectorField> gradi(N_Ori);

for(label count = 0; count < N_Ori; count++)
{

    word name_gradi ("gradi." + name(count));

    gradi.set
    (
        count,
        new volVectorField
        (
            IOobject
            (
                name_gradi,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,//IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("gradi",dimensionSet(0, 0, 0, 0, 0),vector::zero), // NEED TO SET DIMENSIONS TO CHEM POT AFTER I DEBUG
            // zeroGradientFvPatchVectorField::typeName
             fixedValueFvPatchVectorField::typeName
        )
    );
}


tensor rot[N_Ori]={tensor::zero};

for (label i = 0; i < N_Ori; i++)
{


 tensor temprot(
    Foam::cos(phi2[i])*Foam::cos(phi3[i]), 
    (Foam::sin(phi1[i])*Foam::sin(phi2[i])*Foam::cos(phi3[i]))-(Foam::cos(phi1[i])*Foam::sin(phi3[i])),//-Foam::cos(phi3)*Foam::sin(phi2), 
    (Foam::cos(phi1[i])*Foam::sin(phi2[i])*Foam::cos(phi3[i]))+(Foam::sin(phi1[i])*Foam::sin(phi3[i])),//Foam::sin(phi2)*Foam::sin(phi3),
   
    Foam::cos(phi2[i])*Foam::sin(phi3[i]), 
    (Foam::sin(phi1[i])*Foam::sin(phi2[i])*Foam::sin(phi3[i]))+(Foam::cos(phi1[i])*Foam::cos(phi3[i])), 
    (Foam::cos(phi1[i])*Foam::sin(phi2[i])*Foam::sin(phi3[i]))-(Foam::sin(phi1[i])*Foam::cos(phi3[i])),
    
    -Foam::sin(phi2[i]), 
    Foam::sin(phi1[i])*Foam::cos(phi2[i]), 
    Foam::cos(phi1[i])*Foam::cos(phi2[i])
);//https://en.wikipedia.org/wiki/Rotation_matrix - 
//"represents an extrinsic rotation whose (improper) Euler angles are α, β, γ, about axes x, y, z" - phi1 ==alpha, phi2==beta,phi3==gamma

rot[i]=temprot;



}



// Need to check if all procs have the same seed for the random gen
scalar randGenSeed = clock::getTime();
reduce(randGenSeed,minOp<scalar>());
Random randObjq(randGenSeed);

quaternion rot2[N_Ori]={quaternion::zero};

for (label i = 0; i < N_Ori; i++)
{

    scalar X0 = randObjq.sample01<scalar>();
scalar X1 = randObjq.sample01<scalar>();
scalar X2 = randObjq.sample01<scalar>();

    scalar theta1 = 2*M_PI*X1;
    scalar theta2 = 2*M_PI*X2;
    scalar s1 = Foam::sin(theta1);
    scalar c1 = Foam::cos(theta1);
    scalar s2 = Foam::sin(theta2);
    scalar c2 = Foam::cos(theta2);
    scalar r1 = Foam::sqrt(1-X0);
    scalar r2 = Foam::sqrt(X0);

    
    // q.normalize();

    // if(mesh.nGeometricD() == 2){
    // quaternion q(s1*r1,vector(c1*r1,s2*r2,c2*r2));
    // rot2[i]=q;
    // }

    // else{
    quaternion q(s1*r1,vector(c1*r1,s2*r2,c2*r2));
    rot2[i]=q;
    // }

    // Info<<q<<endl;

}

quaternion qZero = {quaternion::zero};

Info << "Creating orientation fields\n\tScalar component qw" << endl;
volScalarField qw 
(
    IOobject
    (
        "qw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("qw", dimless, qZero.w())
);

// Info << "\tVector component qv" << endl;
volVectorField qv
(
    IOobject
    (
        "qv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("qv", dimless, qZero.v())
);

// Info<<rot2[0]<<endl;
// Info<<rot2[1]<<endl;
// Info<<rot2[2]<<endl;

// exit(0);


PtrList<volScalarField> phi_Mob(N_Ori);


    for(label count = 0; count < N_Ori; count++)
{

    word name_ni_mob ("n." + name(count)+ "mob");

    phi_Mob.set
    (
        count,
        new volScalarField 
    (
            IOobject
            (
                name_ni_mob,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,//IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh,
	    dimensionedScalar("phi_Mob",dimensionSet(0, 0, 0, 0, 0),0.0),
	    zeroGradientFvPatchScalarField::typeName
        )
    );
}

volScalarField numerator
(
    IOobject
    (
        "numerator",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("numerator",dimensionSet(0,0,0,0,0),0.0)
);


volScalarField denominator
(
    IOobject
    (
        "denominator",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("denominator",dimensionSet(0,0,0,0,0),0.0)
);

volScalarField Aniso_Mob
(
    IOobject
    (
        "Aniso_Mob",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Aniso_Mob",dimensionSet(0,0,0,0,0),0.0)
);

dimensionedScalar lengthcorr ("lengthcorr", dimensionSet(0,1,0,0,0,0,0),1.0);
dimensionedScalar deltaN2=lengthcorr*(1e-8/pow(average(mesh.V()), 1.0/3.0));



dimensionedScalar Aniso_strength("Aniso_strength", dimensionSet(0, 0, 0, 0, 0), readScalar(PhaseFieldProperties.lookup("Aniso_strength")));

volScalarField Tu
(
    IOobject
    (
        "Tu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    T - TLiquidus
);

//microstructure stuff












typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    > transportModelType;

autoPtr<transportModelType> turbulence
(
    new transportModelType
    (
        rho, U, phi, rhoPhi, mixture
    )
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

// Create laser heat source
laserHeatSource laser(mesh);

#include "updateProps.H"



volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh*rhok
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh*rhok;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"


#include "isoAdvector/createPorosity.H"
isoAdvection advector(alpha1, phi, U);
