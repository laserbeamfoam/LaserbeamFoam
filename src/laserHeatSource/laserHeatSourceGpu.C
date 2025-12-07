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

#ifdef FOAM_USE_CUDA
    #include "cuda_runtime.h"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * *  Helper Structs * * * * * * * * * * * * * * * * //

#ifdef FOAM_USE_CUDA

// Simple POD representation of a ray for device transfers
struct DeviceRay
{
    double3 pos;
    double3 dir;
    double  power;
    int     currentCell;
    int     globalIndex;
};

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


// -------------------------------------------------------------------------
// Stub GPU implementation: just call the CPU backend
// -------------------------------------------------------------------------

void laserHeatSource::updateDepositionGPU
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
    // --- 1) Access the mesh and cell volumes ---
    const fvMesh& mesh = deposition_.mesh();
    const scalarField& VI = mesh.V();
    const label nCells = VI.size();

    Info<< "laserHeatSource::updateDepositionGPU: nCells = "
        << nCells << endl;

    if (nCells > 0)
    {
        // Host-side buffer to receive data back
        List<scalar> VIcopy(nCells, scalar(0));

        // Device pointer
        scalar* dVI = nullptr;

        // --- 2) Allocate on device ---
        cudaError_t err =
            cudaMalloc(reinterpret_cast<void**>(&dVI),
                       nCells*sizeof(scalar));

        if (err != cudaSuccess)
        {
            Info<< "laserHeatSource::updateDepositionGPU: cudaMalloc(VI) failed: "
                << cudaGetErrorString(err)
                << ". Falling back to CPU backend." << endl;
        }
        else
        {
            // --- 3) Copy VI -> device ---
            err = cudaMemcpy
            (
                dVI,
                VI.begin(),
                nCells*sizeof(scalar),
                cudaMemcpyHostToDevice
            );

            if (err != cudaSuccess)
            {
                Info<< "laserHeatSource::updateDepositionGPU: cudaMemcpy H2D(VI) failed: "
                    << cudaGetErrorString(err)
                    << endl;
            }
            else
            {
                // --- 4) Copy VI back to host (device -> VIcopy) ---
                err = cudaMemcpy
                (
                    VIcopy.begin(),
                    dVI,
                    nCells*sizeof(scalar),
                    cudaMemcpyDeviceToHost
                );

                if (err != cudaSuccess)
                {
                    Info<< "laserHeatSource::updateDepositionGPU: cudaMemcpy D2H(VI) failed: "
                        << cudaGetErrorString(err)
                        << endl;
                }
                else
                {
                    // Simple sanity check: compare the first element
                    if (nCells > 0)
                    {
                        Info<< "laserHeatSource::updateDepositionGPU: VI[0] = "
                            << VI[0]
                            << ", VIcopy[0] = " << VIcopy[0]
                            << endl;
                    }
                }
            }

            // --- 5) Free device memory ---
            cudaError_t freeErr = cudaFree(dVI);
            if (freeErr != cudaSuccess)
            {
                Info<< "laserHeatSource::updateDepositionGPU: cudaFree(VI) failed: "
                    << cudaGetErrorString(freeErr)
                    << endl;
            }
        }
    }
    else
    {
        Info<< "laserHeatSource::updateDepositionGPU: nCells = 0, "
            << "skipping VI copy test." << endl;
    }


    // ------------------------------
    // 2) Ray creation + copy test
    // ------------------------------
    {
        Info<< "laserHeatSource::updateDepositionGPU: creating test rays "
            << "for device transfer." << endl;

        // Reuse existing ray initialisation
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
            currentLaserPower,  // Q_cond value
            laserRadius         // beam_radius (a_cond.value())
        );

        const label nRays = rays.size();
        Info<< "laserHeatSource::updateDepositionGPU: nRays = "
            << nRays << endl;

        if (nRays > 0)
        {
            // Host-side flattened representation
            List<DeviceRay> hRays(nRays);

            for (label i = 0; i < nRays; ++i)
            {
                const compactRay& r = rays[i];
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
                    << ". Skipping ray copy test." << endl;
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
                    }
                }

                cudaError_t freeErr = cudaFree(dRays);
                if (freeErr != cudaSuccess)
                {
                    Info<< "laserHeatSource::updateDepositionGPU: cudaFree(dRays) failed: "
                        << cudaGetErrorString(freeErr) << endl;
                }
            }
        }
        else
        {
            Info<< "laserHeatSource::updateDepositionGPU: no rays created; "
                << "skipping ray copy test." << endl;
        }
    }
#else
    Info<< "laserHeatSource::updateDepositionGPU: compiled without CUDA "
        << "support -> using CPU backend." << endl;
#endif

    // For now, always call the CPU implementation
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
}

} // End namespace Foam
