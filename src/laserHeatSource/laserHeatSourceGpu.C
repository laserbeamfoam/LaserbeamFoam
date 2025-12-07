/*--------------------------------*- C++ -*----------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "laserHeatSource.H"
#include "fvc.H"
#include "constants.H"

#ifdef FOAM_USE_CUDA
    #include "cuda_runtime.h"
    #include <thrust/device_vector.h>
    #include <thrust/device_ptr.h>
    #include <thrust/copy.h>
    #include "deviceRayConversion.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * *  Helper Structs * * * * * * * * * * * * * * * * //

#ifdef FOAM_USE_CUDA


// C-linkage function implemented in liblaserHeatSourceCuda.so
extern "C"
void launchGpuTraceRayKernel
(
    DeviceRay*       dRays,
    int              nRays,
    const double*    dVI,
    const double*    dAlpha,
    const double3*   dN,
    const double*    dResistivity,
    int              nCells,
    double3          bbMin,
    double3          bbMax,
    int3             voxelDims,
    double3          voxelSize,
    const int*       dVoxelCell,
    const double3*   dCellMin,
    const double3*   dCellMax,
    double           plasmaFreq,
    double           omega,
    double           dep_cutoff,
    double           rayPowerAbsTol,
    double*          dDeposition,
    int              maxSteps
);


#endif


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //


bool laserHeatSource::detectGPU()
{
#ifdef FOAM_USE_CUDA
    Info<< "laserHeatSource::detectGPU: compiled with CUDA support "
        << "(FOAM_USE_CUDA is defined)." << endl;

    int nDevices = 0;
    cudaError_t err = cudaGetDeviceCount(&nDevices);

    if (err != cudaSuccess)
    {
        Info<< "laserHeatSource::detectGPU: cudaGetDeviceCount failed: "
            << cudaGetErrorString(err)
            << " -> disabling GPU backend." << endl;
        return false;
    }

    Info<< "laserHeatSource::detectGPU: cudaGetDeviceCount reports "
        << nDevices << " device(s)." << endl;

    if (nDevices <= 0)
    {
        Info<< "laserHeatSource::detectGPU: no CUDA devices found -> "
            << "disabling GPU backend." << endl;
        return false;
    }

    const int devId = Pstream::myProcNo() % nDevices;
    err = cudaSetDevice(devId);

    if (err != cudaSuccess)
    {
        Info<< "laserHeatSource::detectGPU: cudaSetDevice(" << devId
            << ") failed: " << cudaGetErrorString(err)
            << " -> disabling GPU backend." << endl;
        return false;
    }

    Info<< "laserHeatSource::detectGPU: using CUDA device " << devId
        << " on rank " << Pstream::myProcNo() << endl;

    return true;
#else
    Info<< "laserHeatSource::detectGPU: compiled WITHOUT CUDA support "
        << "(FOAM_USE_CUDA is NOT defined) -> GPU disabled." << endl;
    return false;
#endif
}


void Foam::laserHeatSource::propagateRaysGPU
(
    DynamicList<compactRay>& remainingGlobalRays,
    const fvMesh& mesh,
    const scalarField& VI,
    const volScalarField& alphaFiltered,
    const volVectorField& nFiltered,
    const volScalarField& resistivity_in,
    const scalar plasma_frequency,
    const scalar angular_frequency,
    const scalar dep_cutoff,
    const scalar rayPowerAbsTol,
    const boundBox& globalBB,
    const label maxLocalSearch, // unused in GPU path
    const label laserID
)
{
#ifndef FOAM_USE_CUDA

    Info<< "propagateRaysGPU: CUDA not enabled (FOAM_USE_CUDA not defined); "
        << "falling back to CPU implementation." << nl;

    propagateRaysCPU
    (
        remainingGlobalRays,
        mesh,
        VI,
        alphaFiltered,
        nFiltered,
        resistivity_in,
        plasma_frequency,
        angular_frequency,
        dep_cutoff,
        rayPowerAbsTol,
        globalBB,
        maxLocalSearch,
        laserID
    );
    return;

#else

    // For now, only support single-rank + voxel mode on GPU
    if (!voxelInitialised_ || cellSearchMode_ == CSM_LEGACY)
    {
        Info<< "propagateRaysGPU: voxel data not initialised or legacy "
            << "cell search in use; falling back to CPU." << nl;

        propagateRaysCPU
        (
            remainingGlobalRays,
            mesh,
            VI,
            alphaFiltered,
            nFiltered,
            resistivity_in,
            plasma_frequency,
            angular_frequency,
            dep_cutoff,
            rayPowerAbsTol,
            globalBB,
            maxLocalSearch,
            laserID
        );
        return;
    }

    if (Pstream::nProcs() > 1)
    {
        Info<< "propagateRaysGPU: MPI run detected (nProcs = "
            << Pstream::nProcs()
            << "): GPU ray transfer across ranks not implemented yet; "
            << "falling back to CPU." << nl;

        propagateRaysCPU
        (
            remainingGlobalRays,
            mesh,
            VI,
            alphaFiltered,
            nFiltered,
            resistivity_in,
            plasma_frequency,
            angular_frequency,
            dep_cutoff,
            rayPowerAbsTol,
            globalBB,
            maxLocalSearch,
            laserID
        );
        return;
    }

    const label nRays = remainingGlobalRays.size();
    if (nRays == 0)
    {
        Info<< "propagateRaysGPU: no rays to trace; returning." << nl;
        return;
    }

    Info<< "propagateRaysGPU: tracing " << nRays
        << " rays on GPU (single kernel)." << nl;

    const label nCells = VI.size();

    // --------------------------
    // 1) Host-side packing
    // --------------------------

    // Rays -> DeviceRay array
    List<DeviceRay> hRays(nRays);
    forAll(remainingGlobalRays, i)
    {
        hRays[i] = toDeviceRay(remainingGlobalRays[i]);
    }

    // Cell fields
    const scalarField& alphaI = alphaFiltered.internalField();
    const vectorField& nI     = nFiltered.internalField();
    const scalarField& resI   = resistivity_in.internalField();

    // Pack normals into double3
    List<double3> hN(nCells);
    forAll(nI, cellI)
    {
        const vector& nv = nI[cellI];
        hN[cellI] = make_double3(nv.x(), nv.y(), nv.z());
    }

    // Pack cell AABBs into double3
    const label nCellsLocal = cellMin_.size();
    List<double3> hCellMin(nCellsLocal);
    List<double3> hCellMax(nCellsLocal);
    forAll(cellMin_, cellI)
    {
        const vector& cMin = cellMin_[cellI];
        const vector& cMax = cellMax_[cellI];
        hCellMin[cellI] = make_double3(cMin.x(), cMin.y(), cMin.z());
        hCellMax[cellI] = make_double3(cMax.x(), cMax.y(), cMax.z());
    }

    // Deposition accumulator on host (from GPU)
    scalarField gpuDeposition(nCells, 0.0);

    // Voxel info
    const int3 voxelDims = make_int3
    (
        static_cast<int>(nVoxelsX_),
        static_cast<int>(nVoxelsY_),
        static_cast<int>(nVoxelsZ_)
    );

    // Compute voxel sizes from edges (safer than relying on stored dx_/dy_/dz_)
    const scalar dx =
        (xEdges_.last() - xEdges_.first())/scalar(nVoxelsX_);
    const scalar dy =
        (yEdges_.last() - yEdges_.first())/scalar(nVoxelsY_);
    const scalar dz =
        (zEdges_.last() - zEdges_.first())/scalar(nVoxelsZ_);

    const double3 voxelSize = make_double3(dx, dy, dz);

    // --------------------------
    // 2) Device allocations
    // --------------------------

    DeviceRay* dRays = nullptr;
    double* dVI = nullptr;
    double* dAlpha = nullptr;
    double3* dN = nullptr;
    double* dRes = nullptr;
    int* dVoxelCell = nullptr;
    double3* dCellMin = nullptr;
    double3* dCellMax = nullptr;
    double* dDeposition = nullptr;

    auto checkCuda = [&](cudaError_t err, const char* what) -> bool
    {
        if (err != cudaSuccess)
        {
            Info<< "propagateRaysGPU: CUDA error in " << what << ": "
                << cudaGetErrorString(err) << nl;
            return false;
        }
        return true;
    };

    // Rays
    if
    (
        !checkCuda
        (
            cudaMalloc(reinterpret_cast<void**>(&dRays),
                       nRays*sizeof(DeviceRay)),
            "cudaMalloc(dRays)"
        )
    )
    {
        return;
    }

    if
    (
        !checkCuda
        (
            cudaMemcpy
            (
                dRays,
                hRays.begin(),
                nRays*sizeof(DeviceRay),
                cudaMemcpyHostToDevice
            ),
            "cudaMemcpy(dRays)"
        )
    )
    {
        cudaFree(dRays);
        return;
    }

    // VI
    if
    (
        !checkCuda
        (
            cudaMalloc(reinterpret_cast<void**>(&dVI),
                       nCells*sizeof(double)),
            "cudaMalloc(dVI)"
        )
    )
    {
        cudaFree(dRays);
        return;
    }

    if
    (
        !checkCuda
        (
            cudaMemcpy
            (
                dVI,
                VI.begin(),
                nCells*sizeof(double),
                cudaMemcpyHostToDevice
            ),
            "cudaMemcpy(dVI)"
        )
    )
    {
        cudaFree(dRays);
        cudaFree(dVI);
        return;
    }

    // alpha
    if
    (
        !checkCuda
        (
            cudaMalloc(reinterpret_cast<void**>(&dAlpha),
                       nCells*sizeof(double)),
            "cudaMalloc(dAlpha)"
        )
    )
    {
        cudaFree(dRays);
        cudaFree(dVI);
        return;
    }

    if
    (
        !checkCuda
        (
            cudaMemcpy
            (
                dAlpha,
                alphaI.begin(),
                nCells*sizeof(double),
                cudaMemcpyHostToDevice
            ),
            "cudaMemcpy(dAlpha)"
        )
    )
    {
        cudaFree(dRays);
        cudaFree(dVI);
        cudaFree(dAlpha);
        return;
    }

    // normals
    if
    (
        !checkCuda
        (
            cudaMalloc(reinterpret_cast<void**>(&dN),
                       nCells*sizeof(double3)),
            "cudaMalloc(dN)"
        )
    )
    {
        cudaFree(dRays);
        cudaFree(dVI);
        cudaFree(dAlpha);
        return;
    }

    if
    (
        !checkCuda
        (
            cudaMemcpy
            (
                dN,
                hN.begin(),
                nCells*sizeof(double3),
                cudaMemcpyHostToDevice
            ),
            "cudaMemcpy(dN)"
        )
    )
    {
        cudaFree(dRays);
        cudaFree(dVI);
        cudaFree(dAlpha);
        cudaFree(dN);
        return;
    }

    // resistivity
    if
    (
        !checkCuda
        (
            cudaMalloc(reinterpret_cast<void**>(&dRes),
                       nCells*sizeof(double)),
            "cudaMalloc(dRes)"
        )
    )
    {
        cudaFree(dRays);
        cudaFree(dVI);
        cudaFree(dAlpha);
        cudaFree(dN);
        return;
    }

    if
    (
        !checkCuda
        (
            cudaMemcpy
            (
                dRes,
                resI.begin(),
                nCells*sizeof(double),
                cudaMemcpyHostToDevice
            ),
            "cudaMemcpy(dRes)"
        )
    )
    {
        cudaFree(dRays);
        cudaFree(dVI);
        cudaFree(dAlpha);
        cudaFree(dN);
        cudaFree(dRes);
        return;
    }

    // voxelCell
    const label nVoxels = voxelCell_.size();
    if
    (
        !checkCuda
        (
            cudaMalloc(reinterpret_cast<void**>(&dVoxelCell),
                       nVoxels*sizeof(int)),
            "cudaMalloc(dVoxelCell)"
        )
    )
    {
        cudaFree(dRays);
        cudaFree(dVI);
        cudaFree(dAlpha);
        cudaFree(dN);
        cudaFree(dRes);
        return;
    }

    if
    (
        !checkCuda
        (
            cudaMemcpy
            (
                dVoxelCell,
                voxelCell_.begin(),
                nVoxels*sizeof(int),
                cudaMemcpyHostToDevice
            ),
            "cudaMemcpy(dVoxelCell)"
        )
    )
    {
        cudaFree(dRays);
        cudaFree(dVI);
        cudaFree(dAlpha);
        cudaFree(dN);
        cudaFree(dRes);
        cudaFree(dVoxelCell);
        return;
    }

    // cellMin / cellMax
    if (!checkCuda
        (
            cudaMalloc(reinterpret_cast<void**>(&dCellMin),
                       nCellsLocal*sizeof(double3)),
            "cudaMalloc(dCellMin)"
        ))
    {
        cudaFree(dRays);
        cudaFree(dVI);
        cudaFree(dAlpha);
        cudaFree(dN);
        cudaFree(dRes);
        cudaFree(dVoxelCell);
        return;
    }

    if (!checkCuda
        (
            cudaMalloc(reinterpret_cast<void**>(&dCellMax),
                       nCellsLocal*sizeof(double3)),
            "cudaMalloc(dCellMax)"
        ))
    {
        cudaFree(dRays);
        cudaFree(dVI);
        cudaFree(dAlpha);
        cudaFree(dN);
        cudaFree(dRes);
        cudaFree(dVoxelCell);
        cudaFree(dCellMin);
        return;
    }

    if (!checkCuda
        (
            cudaMemcpy
            (
                dCellMin,
                hCellMin.begin(),
                nCellsLocal*sizeof(double3),
                cudaMemcpyHostToDevice
            ),
            "cudaMemcpy(dCellMin)"
        ))
    {
        cudaFree(dRays);
        cudaFree(dVI);
        cudaFree(dAlpha);
        cudaFree(dN);
        cudaFree(dRes);
        cudaFree(dVoxelCell);
        cudaFree(dCellMin);
        cudaFree(dCellMax);
        return;
    }

    if (!checkCuda
        (
            cudaMemcpy
            (
                dCellMax,
                hCellMax.begin(),
                nCellsLocal*sizeof(double3),
                cudaMemcpyHostToDevice
            ),
            "cudaMemcpy(dCellMax)"
        ))
    {
        cudaFree(dRays);
        cudaFree(dVI);
        cudaFree(dAlpha);
        cudaFree(dN);
        cudaFree(dRes);
        cudaFree(dVoxelCell);
        cudaFree(dCellMin);
        cudaFree(dCellMax);
        return;
    }

    // deposition
    if (!checkCuda
        (
            cudaMalloc(reinterpret_cast<void**>(&dDeposition),
                       nCells*sizeof(double)),
            "cudaMalloc(dDeposition)"
        ))
    {
        cudaFree(dRays);
        cudaFree(dVI);
        cudaFree(dAlpha);
        cudaFree(dN);
        cudaFree(dRes);
        cudaFree(dVoxelCell);
        cudaFree(dCellMin);
        cudaFree(dCellMax);
        return;
    }

    cudaMemset(dDeposition, 0, nCells*sizeof(double));

    // --------------------------
    // 3) Launch kernel
    // --------------------------

    const point& bbMinP = globalBB.min();
    const point& bbMaxP = globalBB.max();

    const double3 bbMin = make_double3(bbMinP.x(), bbMinP.y(), bbMinP.z());
    const double3 bbMax = make_double3(bbMaxP.x(), bbMaxP.y(), bbMaxP.z());

    const int maxSteps = 10000; // safety cap; adjust later if needed

    Info<< "propagateRaysGPU: launching gpuTraceRayKernel (nRays = "
        << nRays << ", maxSteps = " << maxSteps << ")." << nl;

    launchGpuTraceRayKernel
    (
        dRays,
        static_cast<int>(nRays),
        dVI,
        dAlpha,
        dN,
        dRes,
        static_cast<int>(nCells),
        bbMin,
        bbMax,
        voxelDims,
        voxelSize,
        dVoxelCell,
        dCellMin,
        dCellMax,
        plasma_frequency,
        angular_frequency,
        dep_cutoff,
        rayPowerAbsTol,
        dDeposition,
        maxSteps
    );

    // --------------------------
    // 4) Copy results back
    // --------------------------

    // Deposition
    if (!checkCuda
        (
            cudaMemcpy
            (
                gpuDeposition.begin(),
                dDeposition,
                nCells*sizeof(double),
                cudaMemcpyDeviceToHost
            ),
            "cudaMemcpy(dDeposition -> host)"
        ))
    {
        // Even if this fails, free memory below
    }
    else
    {
        forAll(gpuDeposition, cellI)
        {
            deposition_[cellI] += gpuDeposition[cellI];
        }
    }

    // Rays (optional – may all be dead)
    if (checkCuda
        (
            cudaMemcpy
            (
                hRays.begin(),
                dRays,
                nRays*sizeof(DeviceRay),
                cudaMemcpyDeviceToHost
            ),
            "cudaMemcpy(dRays -> host)"
        ))
    {
        remainingGlobalRays.clear();
        forAll(hRays, i)
        {
            const DeviceRay& dR = hRays[i];

            // Keep only rays still alive and inside globalBB
            if
            (
                dR.power > rayPowerAbsTol
             && globalBB.contains(point(dR.pos.x, dR.pos.y, dR.pos.z))
            )
            {
                remainingGlobalRays.append(toCompactRay(dR));
            }
        }
    }

    // --------------------------
    // 5) Free device memory
    // --------------------------
    cudaFree(dRays);
    cudaFree(dVI);
    cudaFree(dAlpha);
    cudaFree(dN);
    cudaFree(dRes);
    cudaFree(dVoxelCell);
    cudaFree(dCellMin);
    cudaFree(dCellMax);
    cudaFree(dDeposition);

    // Note: rayPaths_ are not updated here yet, since we do not track
    // full paths on the GPU. That can be added later if needed.
#endif
}


void Foam::laserHeatSource::updateDepositionGPU
(
    const volScalarField& alphaFiltered,
    const volVectorField& nFiltered,
    const volScalarField& resistivity_in,
    const label laserID,
    const vector& currentLaserPosition,
    const scalar currentLaserPower,
    const scalar laserRadius,
    const label N_sub_divisions,
    const label nRadial,
    const label nAngular,
    const vector& V_incident,
    const scalar wavelength,
    const scalar e_num_density,
    const scalar dep_cutoff,
    const scalar Radius_Flavour,
    const Switch useLocalSearch,
    const label maxLocalSearch,
    const scalar rayPowerRelTol,
    const boundBox& globalBB
)
{
#ifdef FOAM_USE_CUDA
    Info<< "laserHeatSource::updateDepositionGPU: using GPU backend"
        << endl;

    const fvMesh& mesh  = deposition_.mesh();
    //const Time& runTime = mesh.time();
    const scalarField VI = mesh.V();
    const scalar pi = constant::mathematical::pi;

    const dimensionedScalar a_cond
    (
        "a_cond", dimensionSet(0, 1, 0, 0, 0), laserRadius
    );
    const dimensionedScalar Q_cond
    (
        "Q_cond", dimensionSet(1, 2, -3, 0, 0), currentLaserPower
    );

    const scalar plasma_frequency = Foam::sqrt
    (
        (
            e_num_density
           *constant::electromagnetic::e.value()
           *constant::electromagnetic::e.value()
        )
       /(
           constant::atomic::me.value()
          *constant::electromagnetic::epsilon0.value()
       )
    );
    const scalar angular_frequency =
        2.0*pi*constant::universal::c.value()/wavelength;

    if (debug)
    {
        Info<< "useLocalSearch: " << useLocalSearch << nl << nl << endl;
    }

    const scalar beam_radius = a_cond.value();

    // Take references for brevity
    const volVectorField& nFilteredV = nFiltered;
    const volScalarField& alphaFilteredV = alphaFiltered;

    // Create the initial rays (same as CPU path)
    List<compactRay> rays;
    createInitialRays
    (
        rays,
        mesh,
        currentLaserPosition,
        laserRadius,
        N_sub_divisions,
        nRadial,
        nAngular,
        V_incident,
        Radius_Flavour,
        Q_cond.value(),
        beam_radius
    );

    // remainingGlobalRays will store the rays that have yet to propagate
    DynamicList<compactRay> remainingGlobalRays(rays);

    // Reset the ray paths list
    if (Pstream::master())
    {
        rayPaths_[laserID].clear();
        rayPaths_[laserID].setSize(rays.size());
    }

    // Calculate the ray power tolerance as a fraction of the max ray power
    scalar rayPowerAbsTol = 0;
    {
        scalar maxRayPower = 0.0;
        forAll(remainingGlobalRays, rayI)
        {
            maxRayPower = max(maxRayPower, remainingGlobalRays[rayI].power_);
        }

        rayPowerAbsTol = rayPowerRelTol*maxRayPower;

        Info<< "    Max ray power = " << maxRayPower << nl
            << "    Ray power relative tolerance = " << rayPowerRelTol << nl
            << "    Ray power absolute tolerance = " << rayPowerAbsTol << endl;
    }

    Info<< "    Number of rays: "<< remainingGlobalRays.size() << endl;

    // Propagation entry point: currently calls CPU implementation,
    // but separated so we can later switch to a real GPU kernel.
    propagateRaysGPU
    (
        remainingGlobalRays,
        mesh,
        VI,
        alphaFilteredV,
        nFilteredV,
        resistivity_in,
        plasma_frequency,
        angular_frequency,
        dep_cutoff,
        rayPowerAbsTol,
        globalBB,
        maxLocalSearch,
        laserID
    );

    deposition_.correctBoundaryConditions();

    const scalar TotalQ = fvc::domainIntegrate(deposition_).value();
    Info<< "    Total Q deposited (GPU backend): " << TotalQ << endl;

#else
    Info<< "laserHeatSource::updateDepositionGPU: compiled without CUDA "
        << "support -> using pure CPU implementation." << endl;

    // Fallback: just use the CPU path
    updateDepositionCPU
    (
        alphaFiltered,
        nFiltered,
        resistivity_in,
        laserID,
        currentLaserPosition,
        currentLaserPower,
        laserRadius,
        N_sub_divisions,
        nRadial,
        nAngular,
        V_incident,
        wavelength,
        e_num_density,
        dep_cutoff,
        Radius_Flavour,
        useLocalSearch,
        maxLocalSearch,
        rayPowerRelTol,
        globalBB
    );
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
