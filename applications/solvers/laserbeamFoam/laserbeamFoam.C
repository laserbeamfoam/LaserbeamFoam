/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    laserbeamFoam

Group
    grpMultiphaseSolvers

Description
    Ray-Tracing heat source implementation with two phase incompressible VoF
    description of the metallic substrate and shielding gas phase,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.
Authors
    
    Tom Flint, UoM.
    Philip Cardiff, UCD.
    Gowthaman Parivendhan, UCD.
    Joe Robson, UoM.
    Petar Cosic, UCD
    Simon Rodriguez, UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "isoAdvection.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "incompressibleInterPhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "fvcSmooth.H"
#include "dynamicRefineFvMesh.H"

#include "Polynomial.H"
#include "laserHeatSource.H"

// #include "basicKinematicCollidingCloud.H"
// #include "basicKinematicCloud.H"

// For particles with temperature tracking
#include "basicKinematicCloud.H"
#include "basicKinematicParcel.H"

// Add these for particle temperature storage
#include "IOdictionary.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for two incompressible, isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing.\n"
        "With optional mesh motion and mesh topology changes including"
        " adaptive re-meshing."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "MULES/createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    if (interfaceTrackingScheme == "MULES")
    {
        if (!LTS)
        {
            #include "MULES/CourantNo.H"
            #include "setInitialDeltaT.H"
        }
    }
    else if (interfaceTrackingScheme == "isoAdvector")
    {
        #include "isoAdvector/porousCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {





        #include "readControls.H"


    #include "readDyMControls.H"
    






































        

        if (interfaceTrackingScheme == "MULES")
        {
            if (LTS)
            {
                #include "MULES/setRDeltaT.H"
            }
            else
            {
                #include "MULES/CourantNo.H"
                #include "MULES/alphaCourantNo.H"
                #include "MULES/setDeltaT.H"
            }
        } 
        else if (interfaceTrackingScheme == "isoAdvector")
        {
            #include "isoAdvector/porousCourantNo.H"
            #include "isoAdvector/porousAlphaCourantNo.H"
            #include "isoAdvector/setDeltaT.H"
        }
        
        ++runTime;

     

        Info<< "Time = " << runTime.timeName() << nl << endl;
















































// // Reset sources
alphaMeltSource *= 0.0;
particleEnthalpySource *= 0.0;
massSource *= 0.0;

label nMelted = 0;
label nPartiallyMelted = 0;
label nLimited = 0;  // Track how many were limited by space

Info<< "Evolving " << parcels.name() << endl;
parcels.evolve();

// Track how much alpha we're trying to add to each cell
scalarField proposedAlphaAddition(mesh.nCells(), 0.0);
scalarField proposedMassAddition(mesh.nCells(), 0.0);
scalarField proposedEnthalpySink(mesh.nCells(), 0.0);

// First pass: Calculate what each particle wants to add
forAllIter(basicKinematicCloud, parcels, pIter)
{
    basicKinematicParcel& p = pIter();
    
    if (!p.active())
    {
        continue;
    }
    
    label celli = p.cell();
    
    if (celli < 0 || celli >= mesh.nCells())
    {
        continue;
    }
    
    scalar tempP = T[celli];
    
    if (tempP > particle_Tmelt.value())
    {
        scalar originalMass = p.mass();
        label particleID = p.origId();
        label particleProc = p.origProc();
        label particleKey = particleProc * 1000000 + particleID;
        
        scalar remainingFraction = 1.0;
        if (particleRemainingMassFraction.found(particleKey))
        {
            remainingFraction = particleRemainingMassFraction[particleKey];
        }
        
        scalar remainingMass = originalMass * remainingFraction;
        scalar dt = runTime.deltaTValue();
        
        if (dt < SMALL || remainingMass < SMALL)
        {
            if (remainingMass < SMALL)
            {
                p.active(false);
                particleRemainingMassFraction.erase(particleKey);
            }
            continue;
        }
        
        // Calculate desired melting
        scalar meltingTimescale = time_melt.value();
        scalar maxMeltRatePerSecond = originalMass / meltingTimescale;
        scalar maxMeltThisStep = maxMeltRatePerSecond * dt;
        scalar meltedMass = min(remainingMass, maxMeltThisStep);
        

        
        scalar volumeAsMetal = meltedMass / rho1.value();


        scalar cellVolume = mesh.V()[celli];
        scalar deltaAlpha = volumeAsMetal / cellVolume;
        
        // Accumulate proposed additions
        proposedAlphaAddition[celli] += deltaAlpha;
        proposedMassAddition[celli] += meltedMass;
        
        // Energy calculation
        scalar T_inject = Tinject.value();  // ← From dictionary
        scalar Cp_particle = CpParticle.value();  // ← From dictionary
        scalar sensibleHeat = meltedMass * Cp_particle * (tempP - T_inject);
        scalar latentHeat = meltedMass * LfParticle.value();  // ← From dictionary
        
        proposedEnthalpySink[celli] += (sensibleHeat + latentHeat);
    }
}

// Second pass: Apply limits and actually update sources
forAll(proposedAlphaAddition, celli)
{
    if (proposedAlphaAddition[celli] > SMALL)
    {
        scalar currentAlpha = alpha1[celli];
        scalar availableSpace = max(1.0 - currentAlpha, 0.0);
        
        // Limit to 95% of available space for safety
        scalar maxSafeAddition = 0.95 * availableSpace;
        
        scalar actualAlphaAddition = min(proposedAlphaAddition[celli], maxSafeAddition);
        scalar limitFactor = 1.0;
        
        if (proposedAlphaAddition[celli] > maxSafeAddition)
        {
            limitFactor = maxSafeAddition / proposedAlphaAddition[celli];
            nLimited++;
            
            if (runTime.outputTime())
            {
                Info<< "Cell " << celli << " limited: alpha=" << currentAlpha
                    << ", proposed=" << proposedAlphaAddition[celli]
                    << ", actual=" << actualAlphaAddition << endl;
            }
        }
        
        scalar dt = runTime.deltaTValue();
        scalar cellVolume = mesh.V()[celli];
        
        // Apply limited sources
        alphaMeltSource[celli] = actualAlphaAddition / dt;
        massSource[celli] = (proposedMassAddition[celli] * limitFactor) / (cellVolume * dt);
        particleEnthalpySource[celli] = -(proposedEnthalpySink[celli] * limitFactor) / (cellVolume * dt);
    }
}

// Third pass: Update particle states based on what was actually melted
forAllIter(basicKinematicCloud, parcels, pIter)
{
    basicKinematicParcel& p = pIter();
    
    if (!p.active())
    {
        continue;
    }
    
    label celli = p.cell();
    
    if (celli < 0 || celli >= mesh.nCells())
    {
        continue;
    }
    
    scalar tempP = T[celli];
    
    if (tempP > 2000.0)
    {
        scalar originalMass = p.mass();
        label particleID = p.origId();
        label particleProc = p.origProc();
        label particleKey = particleProc * 1000000 + particleID;
        
        scalar remainingFraction = 1.0;
        if (particleRemainingMassFraction.found(particleKey))
        {
            remainingFraction = particleRemainingMassFraction[particleKey];
        }
        
        scalar remainingMass = originalMass * remainingFraction;
        scalar dt = runTime.deltaTValue();
        
        if (dt < SMALL || remainingMass < SMALL)
        {
            continue;
        }
        
        // Calculate what was actually melted (with limiting)
        scalar meltingTimescale = 0.0001;
        scalar maxMeltRatePerSecond = originalMass / meltingTimescale;
        scalar maxMeltThisStep = maxMeltRatePerSecond * dt;
        scalar desiredMelt = min(remainingMass, maxMeltThisStep);
        
        // Apply same limit factor as was used for this cell
        scalar limitFactor = 1.0;
        if (proposedAlphaAddition[celli] > SMALL)
        {
            scalar currentAlpha = alpha1[celli];
            scalar availableSpace = max(1.0 - currentAlpha, 0.0);
            scalar maxSafeAddition = 0.95 * availableSpace;
            
            if (proposedAlphaAddition[celli] > maxSafeAddition)
            {
                limitFactor = maxSafeAddition / proposedAlphaAddition[celli];
            }
        }
        
        scalar actualMeltedMass = desiredMelt * limitFactor;
        scalar newRemainingFraction = (remainingMass - actualMeltedMass) / originalMass;
        
        if (newRemainingFraction < 0.01)
        {
            p.active(false);
            particleRemainingMassFraction.erase(particleKey);
            nMelted++;
        }
        else
        {
            particleRemainingMassFraction.set(particleKey, newRemainingFraction);
            nPartiallyMelted++;
        }
    }
}

// Parallel reductions
reduce(nMelted, sumOp<label>());
reduce(nPartiallyMelted, sumOp<label>());
reduce(nLimited, sumOp<label>());

if (nMelted > 0 || nPartiallyMelted > 0)
{
    scalar totalMassAdded = gSum(massSource.primitiveField() * mesh.V()) 
                          * runTime.deltaTValue();
    
    Info<< "Fully melted: " << nMelted << " particles, "
        << "partially melted: " << nPartiallyMelted << " particles, "
        << "limited by space: " << nLimited << " cells, "
        << "total mass added: " << totalMassAdded*1e3 << " g" << endl;
}

// Check for overfilling
scalar maxAlphaSource = gMax(alphaMeltSource);
scalar maxPredictedAlpha = gMax(alpha1.primitiveField() + alphaMeltSource.primitiveField() * runTime.deltaTValue());

if (maxAlphaSource > SMALL)
{
    Info<< "Max alpha source: " << maxAlphaSource << " 1/s" << endl;
    Info<< "Max predicted alpha after source: " << maxPredictedAlpha << endl;
    
    if (maxPredictedAlpha > 1.0)
    {
        WarningInFunction
            << "Alpha may exceed 1.0! Max predicted = " << maxPredictedAlpha << endl;
    }
}







scalar maxAllowedTempDrop = 200.0;  // Max K drop per timestep
scalar dt = runTime.deltaTValue();

forAll(particleEnthalpySource, celli)
{
    if (particleEnthalpySource[celli] < -SMALL)
    {
        // Maximum cooling power that would drop T by maxAllowedTempDrop
        scalar maxCoolingPower = rhoCp[celli] * maxAllowedTempDrop / dt;
        
        // Limit the source
        if (mag(particleEnthalpySource[celli]) > maxCoolingPower)
        {
            // scalar oldSource = particleEnthalpySource[celli];
            particleEnthalpySource[celli] = -maxCoolingPower;
            
            // // Also need to scale back the mass/alpha sources proportionally
            // // to maintain energy consistency
            // scalar scaleFactor = mag(particleEnthalpySource[celli]) / mag(oldSource);
            // massSource[celli] *= scaleFactor;
            // alphaMeltSource[celli] *= scaleFactor;
        }
    }
}
Info<<"min enthalpy source "<<gMin(particleEnthalpySource)<<endl;
// // fvc::smooth(particleEnthalpySource, 2);






























        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {

            parcels.storeGlobalPositions();

            if (interfaceTrackingScheme == "MULES")
            {
                #include "MULES/firstIter.H"
                #include "MULES/alphaControls.H"
                #include "MULES/alphaEqnSubCycle.H"
            } 
            else if (interfaceTrackingScheme == "isoAdvector")
            {
                #include "isoAdvector/firstIter.H"
                #include "isoAdvector/alphaControls.H"
                #include "isoAdvector/alphaEqnSubCycle.H"
            }

            

            #include "updateProps.H"

            // Update the laser deposition field
            laser.updateDeposition
            (
                alpha_filtered, n_filtered, electrical_resistivity
            );

            mixture.correct();

            if (pimple.frozenFlow())
            {
                continue;
            }

            #include "UEqn.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }


        {// Check the cells that have melted
            volScalarField alphaMetal = 
            mesh.lookupObject<volScalarField>("alpha.metal");
            condition = pos(alphaMetal - 0.5) * pos(epsilon1 - 0.5);
            meltHistory += condition;
        }





        runTime.write();
        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
