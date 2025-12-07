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
void launchGpuRayStepKernel
(
    DeviceRay* dRays,
    int nRays,
    const double* dVI,
    int nCells,
    double bbMinX, double bbMinY, double bbMinZ,
    double bbMaxX, double bbMaxY, double bbMaxZ,
    double rayPowerAbsTol
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
    const label maxLocalSearch,
    const label laserID
)
{
#ifndef FOAM_USE_CUDA
    Info<< "propagateRaysGPU: compiled without CUDA support -> "
        << "using CPU propagation." << endl;

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
    Info<< "propagateRaysGPU: using GPU stepping for iterator distance, "
        << "CPU for cell search and deposition." << endl;

    const volScalarField& alphaFilteredI = alphaFiltered;
    const volVectorField& nFilteredI     = nFiltered;

    // This is a GPU-aware variant of propagateRaysCPU:
    //  - outer while and deposition logic are the same
    //  - the actual "step along ray direction" is done on the GPU.

    while (remainingGlobalRays.size() > 0)
    {
        // Find all rays on the current processor
        DynamicList<compactRay> localRays;

        forAll(remainingGlobalRays, rayI)
        {
            compactRay& curRay = remainingGlobalRays[rayI];

            if
            (
                globalBB.contains(curRay.position_)
             && curRay.power_ > rayPowerAbsTol
            )
            {
                const label myCellID =
                    findCellForRay
                    (
                        curRay.position_,
                        curRay.currentCell_,
                        mesh,
                        maxLocalSearch,
                        debug
                    );

                if (myCellID != -1)
                {
                    // Initialise current cell if needed
                    curRay.currentCell_ = myCellID;
                    localRays.append(curRay);
                }
            }
        }

        if (localRays.size() == 0)
        {
            // No rays left on this processor with sufficient power
            remainingGlobalRays.clear();
            Pstream::combineGather(remainingGlobalRays, combineRayLists());
            Pstream::broadcast(remainingGlobalRays);
            break;
        }

        // Perform one iterator-distance step + culling on the GPU
        applyGpuRayTransform(localRays, VI, globalBB, rayPowerAbsTol);

        // Now do the "physics" part on the CPU, starting from the new positions
        forAll(localRays, rayI)
        {
            compactRay& curRay = localRays[rayI];

            // Find the cell the ray is now in (after GPU step)
            label myCellID =
                findCellForRay
                (
                    curRay.position_,
                    curRay.currentCell_,
                    mesh,
                    maxLocalSearch,
                    debug
                );

            while (myCellID != -1)
            {
                // Update the ray's cellID
                curRay.currentCell_ = myCellID;

                // Update visualisation fields
                rayQ_[myCellID]      += curRay.power_;
                rayNumber_[myCellID]  = curRay.globalRayIndex_;

                if (curRay.power_ < SMALL)
                {
                    curRay.path_.append(curRay.position_);
                    break;
                }
                else if
                (
                    mag(nFilteredI[myCellID]) > 0.5
                 && alphaFilteredI[myCellID] >= dep_cutoff
                )
                {
                    // --- Interface: deposit + reflect (same as CPU) ---

                    const scalar damping_frequency =
                        plasma_frequency*plasma_frequency
                       *constant::electromagnetic::epsilon0.value()
                       *resistivity_in[myCellID];

                    const scalar e_r =
                        1.0
                      - (
                            sqr(plasma_frequency)
                           /(sqr(angular_frequency) + sqr(damping_frequency))
                        );

                    const scalar e_i =
                        (damping_frequency/angular_frequency)
                       *(
                            sqr(plasma_frequency)
                           /(sqr(angular_frequency) + sqr(damping_frequency))
                        );

                    const scalar ref_index =
                        Foam::sqrt
                        (
                            (Foam::sqrt(e_r*e_r + e_i*e_i) + e_r)/2.0
                        );

                    const scalar ext_coefficient =
                        Foam::sqrt
                        (
                            (Foam::sqrt(e_r*e_r + e_i*e_i) - e_r)/2.0
                        );

                    vector n = nFilteredI[myCellID];
                    n /= mag(n);

                    vector d   = curRay.direction_;
                    const scalar dMag = mag(d);

                    if (dMag <= SMALL)
                    {
                        curRay.power_ = 0.0;
                        curRay.path_.append(curRay.position_);
                        break;
                    }

                    d /= dMag;
                    const vector kin = -d;

                    scalar cosTheta = kin & n;

                    if (cosTheta < 0.0)
                    {
                        n        = -n;
                        cosTheta = -cosTheta;
                    }

                    cosTheta =
                        Foam::max
                        (
                            Foam::min(cosTheta, scalar(1.0)),
                            scalar(0.0)
                        );

                    const scalar theta_in = std::acos(cosTheta);
                    const scalar sinTheta = Foam::sin(theta_in);

                    const scalar alpha_laser =
                        Foam::sqrt
                        (
                            Foam::sqrt
                            (
                                sqr
                                (
                                    sqr(ref_index)
                                  - sqr(ext_coefficient)
                                  - sqr(sinTheta)
                                )
                              + 4.0*sqr(ref_index)*sqr(ext_coefficient)
                            )
                          + sqr(ref_index)
                          - sqr(ext_coefficient)
                          - sqr(sinTheta)/2.0
                        );

                    const scalar beta_laser =
                        Foam::sqrt
                        (
                            (
                                Foam::sqrt
                                (
                                    sqr
                                    (
                                        sqr(ref_index)
                                      - sqr(ext_coefficient)
                                      - sqr(sinTheta)
                                    )
                                  + 4.0*sqr(ref_index)*sqr(ext_coefficient)
                                )
                              - sqr(ref_index)
                              + sqr(ext_coefficient)
                              + sqr(sinTheta)
                            )/2.0
                        );

                    const scalar cosTheta_in = Foam::cos(theta_in);

                    scalar R_s =
                        (
                            sqr(alpha_laser)
                          + sqr(beta_laser)
                          - 2.0*alpha_laser*cosTheta_in
                          + sqr(cosTheta_in)
                        )
                       /(
                            sqr(alpha_laser)
                          + sqr(beta_laser)
                          + 2.0*alpha_laser*cosTheta_in
                          + sqr(cosTheta_in)
                        );

                    scalar R_p =
                        R_s
                       *(
                            sqr(alpha_laser)
                          + sqr(beta_laser)
                          - 2.0*alpha_laser*sinTheta*Foam::tan(theta_in)
                          + sqr(sinTheta)*sqr(Foam::tan(theta_in))
                        )
                       /(
                            sqr(alpha_laser)
                          + sqr(beta_laser)
                          + 2.0*alpha_laser*sinTheta*Foam::tan(theta_in)
                          + sqr(sinTheta)*sqr(Foam::tan(theta_in))
                        );

                    R_s =
                        Foam::max
                        (
                            Foam::min(R_s, scalar(1.0)), scalar(0.0)
                        );
                    R_p =
                        Foam::max
                        (
                            Foam::min(R_p, scalar(1.0)), scalar(0.0)
                        );

                    const scalar R      = 0.5*(R_s + R_p);
                    scalar absorptivity = 1.0 - R;
                    absorptivity =
                        Foam::max
                        (
                            Foam::min(absorptivity, scalar(1.0)),
                            scalar(0.0)
                        );

                    if (debug)
                    {
                        Info<< "ray " << curRay.globalRayIndex_
                            << ", cell " << myCellID
                            << ", theta_in = " << theta_in
                            << ", R_s = " << R_s
                            << ", R_p = " << R_p
                            << ", absorptivity = " << absorptivity << endl;
                    }

                    deposition_[myCellID] +=
                        absorptivity*curRay.power_/VI[myCellID];

                    curRay.power_ *= (1.0 - absorptivity);

                    const vector dRef =
                        curRay.direction_ - 2.0*(curRay.direction_ & n)*n;

                    curRay.direction_ = dRef;
                }
                else if (alphaFilteredI[myCellID] >= dep_cutoff)
                {
                    if (debug)
                    {
                        Info<< "Bulk absorption at cell " << myCellID
                            << ", alpha = " << alphaFilteredI[myCellID]
                            << ", |n| = " << mag(nFilteredI[myCellID])
                            << ", power = " << curRay.power_ << endl;
                    }

                    deposition_[myCellID] += curRay.power_/VI[myCellID];
                    curRay.power_ = 0.0;
                    curRay.path_.append(curRay.position_);
                    break;
                }

                curRay.path_.append(curRay.position_);

                // After handling this cell, we break the inner while so that
                // the *next* iterator-distance step (and position update) is
                // again handled on the GPU in the next outer iteration.
                break;
            }
        }

        // Sync remaining rays globally
        remainingGlobalRays = localRays;
        Pstream::combineGather(remainingGlobalRays, combineRayLists());
        Pstream::broadcast(remainingGlobalRays);

        if (Pstream::master())
        {
            forAll(remainingGlobalRays, rI)
            {
                const compactRay& curRay = remainingGlobalRays[rI];
                const label rayID        = curRay.globalRayIndex_;
                rayPaths_[laserID][rayID] = curRay.path_;
            }
        }
    }
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


void Foam::laserHeatSource::applyGpuRayTransform
(
    DynamicList<compactRay>& remainingGlobalRays,
    const scalarField& VI,
    const boundBox& globalBB,
    const scalar rayPowerAbsTol
)
{
#ifdef FOAM_USE_CUDA
    const label nRays  = remainingGlobalRays.size();
    const label nCells = VI.size();

    Info<< "laserHeatSource::applyGpuRayTransform: transferring "
        << nRays << " rays to device and back." << endl;

    if (nRays == 0 || nCells == 0)
    {
        return;
    }

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
        dr.step        = 0.0; // currently unused in kernel

        hRays[i] = dr;
    }

    // Device arrays
    DeviceRay* dRays = nullptr;
    double* dVI      = nullptr;

    cudaError_t err =
        cudaMalloc(reinterpret_cast<void**>(&dRays),
                   nRays*sizeof(DeviceRay));

    if (err != cudaSuccess)
    {
        Info<< "laserHeatSource::applyGpuRayTransform: cudaMalloc(dRays) failed: "
            << cudaGetErrorString(err)
            << ". Skipping GPU transform." << endl;
        return;
    }

    err = cudaMalloc(reinterpret_cast<void**>(&dVI),
                     nCells*sizeof(double));

    if (err != cudaSuccess)
    {
        Info<< "laserHeatSource::applyGpuRayTransform: cudaMalloc(dVI) failed: "
            << cudaGetErrorString(err) << endl;
        cudaFree(dRays);
        return;
    }

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
        Info<< "laserHeatSource::applyGpuRayTransform: cudaMemcpy H2D(rays) failed: "
            << cudaGetErrorString(err) << endl;
        cudaFree(dRays);
        cudaFree(dVI);
        return;
    }

    // VI is a scalarField (double in DP build)
    err = cudaMemcpy
    (
        dVI,
        VI.begin(),
        nCells*sizeof(double),
        cudaMemcpyHostToDevice
    );

    if (err != cudaSuccess)
    {
        Info<< "laserHeatSource::applyGpuRayTransform: cudaMemcpy H2D(VI) failed: "
            << cudaGetErrorString(err) << endl;
        cudaFree(dRays);
        cudaFree(dVI);
        return;
    }

    // Launch kernel
    Info<< "laserHeatSource::applyGpuRayTransform: launching GPU iterator "
        << "step kernel on " << nRays << " rays." << endl;

    const point& bbMin = globalBB.min();
    const point& bbMax = globalBB.max();

    launchGpuRayStepKernel
    (
        dRays,
        nRays,
        dVI,
        nCells,
        bbMin.x(), bbMin.y(), bbMin.z(),
        bbMax.x(), bbMax.y(), bbMax.z(),
        rayPowerAbsTol
    );

    cudaError_t kerr = cudaGetLastError();
    if (kerr != cudaSuccess)
    {
        Info<< "laserHeatSource::applyGpuRayTransform: kernel launch error: "
            << cudaGetErrorString(kerr) << endl;
    }

    // Device -> host
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
        Info<< "laserHeatSource::applyGpuRayTransform: cudaMemcpy D2H(rays) failed: "
            << cudaGetErrorString(err) << endl;
        cudaFree(dRays);
        cudaFree(dVI);
        return;
    }

    // Simple debug for the first ray
    {
        const DeviceRay& r0  = hRays[0];
        const DeviceRay& r0c = hRaysCopy[0];

        Info<< "laserHeatSource::applyGpuRayTransform: ray[0].pos before = ("
            << r0.pos.x << ", " << r0.pos.y << ", " << r0.pos.z
            << "), after = ("
            << r0c.pos.x << ", " << r0c.pos.y << ", " << r0c.pos.z
            << ")" << endl;

        Info<< "laserHeatSource::applyGpuRayTransform: ray[0].power before = "
            << r0.power << ", after = " << r0c.power << endl;
    }

    // Map back into remainingGlobalRays
    const label nRemaining = remainingGlobalRays.size();
    const label nToMap     = min(nRays, nRemaining);

    for (label i = 0; i < nToMap; ++i)
    {
        compactRay& r       = remainingGlobalRays[i];
        const DeviceRay& dr = hRaysCopy[i];

        r.position_.x()      = dr.pos.x;
        r.position_.y()      = dr.pos.y;
        r.position_.z()      = dr.pos.z;

        r.direction_.x()     = dr.dir.x;
        r.direction_.y()     = dr.dir.y;
        r.direction_.z()     = dr.dir.z;

        r.power_             = dr.power;
        r.currentCell_       = dr.currentCell;
        r.globalRayIndex_    = dr.globalIndex;
    }

    Info<< "laserHeatSource::applyGpuRayTransform: mapped "
        << nToMap << " device rays back to remainingGlobalRays." << endl;

    cudaFree(dRays);
    cudaFree(dVI);

#else
    (void)remainingGlobalRays;
    (void)VI;
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
