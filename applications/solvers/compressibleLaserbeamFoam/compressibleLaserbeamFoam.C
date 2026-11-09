/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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
    compressibleLaserbeamFoam

Group
    grpMultiphaseSolvers

Description
    cheeky bad solver for wicked simulation of laser-substrate interactions. A true naughty piece of work

Authors
    Tom Flint UoM
    Roman Bialek UoM
    Philip Cardiff UcD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "multiphaseMixtureThermo.H"
#include "turbulentFluidThermoModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "CorrectPhi.H"
#include "cyclicAMIFvPatch.H"

#include "laserHeatSource.H"
#include "mthdModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for N compressible, non-isothermal immiscible fluids"
        " using VOF phase-fraction based interface capturing."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"

    #include "createDynamicFvMesh.H"

    // #include "createMesh.H"
    // #include "createControl.H"
    // #include "createTimeControls.H"

    #include "initContinuityErrs.H"//new
    #include "createDyMControls.H"//new

    // bool bRestartFirstLoop = false; <- declared within createFields.H
    #include "createFields.H"

    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"



    volScalarField& p = mixture.p();
    volScalarField& T = mixture.T();

    turbulence->validate();

    #include "update.H"


    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {

        #include "readMeltingControls.H"
        #include "readDyMControls.H"

        // #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {

            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // AMR merges cells by volume-averaging every field: right for c_i,
                // rhoCv*T, rhoM*U and rhoM*epsilon1, wrong for T, U and epsilon1.
                // Carry energy, momentum and latent heat through it
                volScalarField rhoCvTAmr("rhoCvT.amr", rhoCv*T);
                volVectorField rhoMUAmr("rhoMU.amr", rhoM*U);
                volScalarField rhoMEps1Amr("rhoMEps1.amr", rhoM*epsilon1);

                scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

                mesh.update();

                if (mesh.changing())
                {
                    Info<< "Execution time for mesh.update() = "
                        << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
                        << " s" << endl;

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi && !bRestartFirstLoop)
                    {
                        // Prevent corrections on the first loop of a restarted
                        // solver.

                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.esi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    #include "remapConserved.H"

                    mixture.correct();

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            // Disable flag: this permits correctPhi on subsequent loops.
            bRestartFirstLoop = false;

            // Semi-PIMPLE Implementation - Roman Bialek
            if (pimple.firstIter())
            {
                rho.oldTime() = rho;   // set once per time step
                rhoM.oldTime() = rhoM; // ledger density, likewise

                vDot = mixture.solve(&mass_dot, &vDotDt, &closureSrc);
                vDot.correctBoundaryConditions();
                mass_dot.correctBoundaryConditions();


                #include "update.H"


            }


            rho = mixture.rho();

            // ledger density for UEqn/pEqn/TEqn, frozen over the outer correctors like rhoPhi
            if (pimple.firstIter())
            {
                #include "updateRhoM.H"
            }



            // Ray-trace expensive: evaluate once per time step.
            if (pimple.firstIter())
            {
                laser.updateDeposition
                (
                    condensateFiltered, n_filtered, electrical_resistivity
                );
            }

            // PIMPLE Section ((U-p)-H-T coupling)

            #include "UEqn.H"
            if (mthd.valid())
            {
                mthd->solve(phi, U);
            }
            #include "TEqn.H"
            #include "volumeClosure.H"

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

        mixture.writeOldTimeValues();
        #include "writeOldTimeStorage.H"

        runTime.write();

        // Write ray paths to VTK files
        if (runTime.outputTime())
        {
            laser.writeRayPathsToVTK();
        }

        runTime.printExecutionTime(Info);
    }

    // Write a VTK series file for easy-opening of the ray files
    laser.writeRayPathVTKSeriesFile();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
