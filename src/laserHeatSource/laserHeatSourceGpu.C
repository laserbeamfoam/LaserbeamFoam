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
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * *  Helper Structs * * * * * * * * * * * * * * * * //

#ifdef FOAM_USE_CUDA

// Simple representation of a ray for device transfers
struct DeviceRay
{
    double3 pos;
    double3 dir;
    double  power;
    int     currentCell;
    int     globalIndex;
    double  step;
};

// C-linkage function implemented in liblaserHeatSourceCuda.so
extern "C"
void launchNoopRayKernel(DeviceRay* dRays, int nRays);

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
    const label maxLocalSearch,
    const label laserID
)
{
#ifdef FOAM_USE_CUDA
    Info<< "propagateRaysGPU: stub implementation -> using CPU propagation."
        << endl;
#else
    Info<< "propagateRaysGPU: compiled without CUDA support -> "
        << "using CPU propagation." << endl;
#endif

    // For now just forward to the CPU implementation
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
    Info<< "laserHeatSource::updateDepositionGPU: using GPU backend stub "
        << "(rays copied to device, propagation still on CPU)." << endl;

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

    // Compute a conservative step length (min cell dimension)
    scalar minCellH = GREAT;
    for (label cellI = 0; cellI < VI.size(); ++cellI)
    {
        scalar h = cbrt(VI[cellI]);
        minCellH = min(minCellH, h);
    }

    // Take a fraction of the min cell size
    scalar stepSize = 0.2 * minCellH;

    Info<< "laserHeatSource::updateDepositionGPU: stepSize = "
        << stepSize << endl;

    // ---------------------------------------------------------------------
    // GPU STEP (stub): copy the *actual* physical rays to/from the device
    // ---------------------------------------------------------------------
    {
        const label nRays = remainingGlobalRays.size();

        if (nRays > 0)
        {
            Info<< "laserHeatSource::updateDepositionGPU: transferring "
                << nRays << " rays to device and back." << endl;

            // Host-side flattened representation
            List<DeviceRay> hRays(nRays);

            for (label i = 0; i < nRays; ++i)
            {
                const compactRay& r = remainingGlobalRays[i];
                DeviceRay dr;

                dr.pos = make_double3
                (
                    r.position_.x(),
                    r.position_.y(),
                    r.position_.z()
                );
                dr.dir = make_double3
                (
                    r.direction_.x(),
                    r.direction_.y(),
                    r.direction_.z()
                );
                dr.power       = r.power_;
                dr.currentCell = r.currentCell_;
                dr.globalIndex = r.globalRayIndex_;
                dr.step = stepSize;

                hRays[i] = dr;
            }

            // Device buffer
            DeviceRay* dRays = nullptr;

            cudaError_t err =
                cudaMalloc(reinterpret_cast<void**>(&dRays),
                           nRays*sizeof(DeviceRay));

            if (err != cudaSuccess)
            {
                Info<< "laserHeatSource::updateDepositionGPU: cudaMalloc(dRays) failed: "
                    << cudaGetErrorString(err)
                    << ". Skipping ray device transfer." << endl;
            }
            else
            {
                // Copy host -> device
                err = cudaMemcpy
                (
                    dRays,
                    hRays.begin(),
                    nRays*sizeof(DeviceRay),
                    cudaMemcpyHostToDevice
                );

                if (err != cudaSuccess)
                {
                    Info<< "laserHeatSource::updateDepositionGPU: cudaMemcpy H2D(rays) failed: "
                        << cudaGetErrorString(err) << endl;
                }
                else
                {
                    // --- NEW: launch trivial CUDA kernel on the device rays ---
                    Info<< "laserHeatSource::updateDepositionGPU: launching noop "
                        << "ray kernel on " << nRays << " rays." << endl;

                    launchNoopRayKernel(dRays, nRays);

                    // Optionally check for errors
                    cudaError_t kerr = cudaGetLastError();
                    if (kerr != cudaSuccess)
                    {
                        Info<< "laserHeatSource::updateDepositionGPU: kernel launch error: "
                            << cudaGetErrorString(kerr) << endl;
                    }

                    // Copy back into a second host buffer
                    List<DeviceRay> hRaysCopy(nRays);
                    err = cudaMemcpy
                    (
                        hRaysCopy.begin(),
                        dRays,
                        nRays*sizeof(DeviceRay),
                        cudaMemcpyDeviceToHost
                    );

                    if (err != cudaSuccess)
                    {
                        Info<< "laserHeatSource::updateDepositionGPU: cudaMemcpy D2H(rays) failed: "
                            << cudaGetErrorString(err) << endl;
                    }
                    else
                    {
                        const DeviceRay& r0  = hRays[0];
                        const DeviceRay& r0c = hRaysCopy[0];

                        Info<< "laserHeatSource::updateDepositionGPU: ray[0].pos = ("
                            << r0.pos.x << ", " << r0.pos.y << ", " << r0.pos.z
                            << "), copy = ("
                            << r0c.pos.x << ", " << r0c.pos.y << ", " << r0c.pos.z
                            << ")" << endl;

                        Info<< "laserHeatSource::updateDepositionGPU: ray[0].power = "
                            << r0.power << ", copy = " << r0c.power << endl;

                        // Write the device-updated rays back into
                        // remainingGlobalRays
                        const label nRemaining = remainingGlobalRays.size();
                        const label nToMap = min(nRays, nRemaining);

                        for (label i = 0; i < nToMap; ++i)
                        {
                            compactRay& r = remainingGlobalRays[i];
                            const DeviceRay& dr = hRaysCopy[i];

                            r.position_.x() = dr.pos.x;
                            r.position_.y() = dr.pos.y;
                            r.position_.z() = dr.pos.z;

                            r.direction_.x() = dr.dir.x;
                            r.direction_.y() = dr.dir.y;
                            r.direction_.z() = dr.dir.z;

                            r.power_       = dr.power;
                            r.currentCell_ = dr.currentCell;
                            r.globalRayIndex_ = dr.globalIndex;
                        }

                        Info<< "laserHeatSource::updateDepositionGPU: mapped "
                            << nToMap << " device rays back to remainingGlobalRays."
                            << endl;
                    }
                }
            }
        }
        else
        {
            Info<< "laserHeatSource::updateDepositionGPU: no rays to transfer."
                << endl;
        }
    }

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
    Info<< "    Total Q deposited (GPU backend, CPU propagation): "
        << TotalQ << endl;

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
