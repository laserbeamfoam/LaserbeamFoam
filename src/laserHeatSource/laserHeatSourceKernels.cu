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

#include <cuda_runtime.h>
#include <math.h>
#include <cstdio>
#include "deviceRay.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// laserHeatSourceKernels.cu
#include <cuda_runtime.h>
#include <math.h>

#include "deviceRay.H"   // Only has DeviceRay + cuda_runtime.h, no OpenFOAM

// =====================================================================
// device helper: voxelFindCellForRay (GPU equivalent of findCellVoxelised)
// =====================================================================
__device__ int voxelFindCellForRay(
    const double3& p,
    int currentCell,
    const double3* cellMin,    // size nCells
    const double3* cellMax,    // size nCells
    int nCells,
    double3 bbMin,
    double3 bbMax,
    int3 voxelDims,            // (Nx, Ny, Nz)
    double3 voxelSize,         // (dx, dy, dz)
    const int* voxelCell       // size Nx*Ny*Nz, voxelCell[vI] = cell or -1
)
{
    // 1) Stay in same cell if still inside AABB
    if (currentCell >= 0 && currentCell < nCells)
    {
        const double3 cMin = cellMin[currentCell];
        const double3 cMax = cellMax[currentCell];

        if
        (
            p.x >= cMin.x && p.x <= cMax.x &&
            p.y >= cMin.y && p.y <= cMax.y &&
            p.z >= cMin.z && p.z <= cMax.z
        )
        {
            return currentCell;
        }
    }

    // 2) Map to voxel index
    const double xRel = p.x - bbMin.x;
    const double yRel = p.y - bbMin.y;
    const double zRel = p.z - bbMin.z;

    int i = static_cast<int>(floor(xRel/voxelSize.x));
    int j = static_cast<int>(floor(yRel/voxelSize.y));
    int k = static_cast<int>(floor(zRel/voxelSize.z));

    if
    (
        i < 0 || i >= voxelDims.x ||
        j < 0 || j >= voxelDims.y ||
        k < 0 || k >= voxelDims.z
    )
    {
        return -1;
    }

    const int vIdx = i + voxelDims.x*(j + voxelDims.y*k);
    const int cellI = voxelCell[vIdx];

    return cellI;
}

// =====================================================================
// global kernel: gpuTraceRayKernel
// =====================================================================
__global__ void gpuTraceRayKernel
(
    DeviceRay*       rays,
    int              nRays,
    const double*    VI,           // cell volumes
    const double*    alpha,        // absorption
    const double3*   n,            // normals
    const double*    resistivity,  // resistivity (unused for now)
    int              nCells,
    double3          bbMin,
    double3          bbMax,
    int3             voxelDims,
    double3          voxelSize,
    const int*       voxelCell,
    const double3*   cellMin,
    const double3*   cellMax,
    double           plasmaFreq,       // unused for now (placeholder)
    double           omega,           // unused for now
    double           dep_cutoff,      // unused for now
    double           rayPowerAbsTol,
    double*          deposition,      // size nCells
    int              maxSteps
)
{
    const int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >= nRays) return;

    DeviceRay& r = rays[idx];

    // simple loop to avoid infinite trace
    int steps = 0;
    while (steps < maxSteps)
    {
        // dead ray? stop
        if (r.power <= rayPowerAbsTol) break;

        // outside global bounding box? stop
        if
        (
            r.pos.x < bbMin.x || r.pos.x > bbMax.x ||
            r.pos.y < bbMin.y || r.pos.y > bbMax.y ||
            r.pos.z < bbMin.z || r.pos.z > bbMax.z
        )
        {
            break;
        }

        // find cell
        const int cellI = voxelFindCellForRay(
            r.pos,
            r.currentCell,
            cellMin,
            cellMax,
            nCells,
            bbMin,
            bbMax,
            voxelDims,
            voxelSize,
            voxelCell
        );

        if (cellI < 0 || cellI >= nCells)
        {
            // outside mesh or unknown voxel
            break;
        }

        r.currentCell = cellI;

        const double VIc    = VI[cellI];
        const double alphac = alpha[cellI];

        // iterator distance like CPU
        const double step = (0.5/M_PI)*cbrt(VIc);

        // ---- bulk absorption only (no Fresnel yet) ----
        if (alphac >= dep_cutoff)
        {
            // CPU bulk branch:
            // deposition_[myCellID] += curRay.power_/VI[myCellID];
            // curRay.power_ = 0.0; break;
            const double volSrc = (VIc > 0.0 ? r.power/VIc : 0.0);

            atomicAdd(&deposition[cellI], volSrc);

            r.power = 0.0;
            break;  // ray is fully absorbed in this cell
        }

        // otherwise: just move the ray, no deposition yet
        r.pos.x += step*r.dir.x;
        r.pos.y += step*r.dir.y;
        r.pos.z += step*r.dir.z;


        ++steps;
    }

    // after maxSteps or termination, the ray is done
}

// =====================================================================
// C-linkage launcher (no OpenFOAM types)
// =====================================================================
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
)
{
    dim3 block(256);
    dim3 grid((nRays + block.x - 1)/block.x);

    gpuTraceRayKernel<<<grid, block>>>(
        dRays,
        nRays,
        dVI,
        dAlpha,
        dN,
        dResistivity,
        nCells,
        bbMin,
        bbMax,
        voxelDims,
        voxelSize,
        dVoxelCell,
        dCellMin,
        dCellMax,
        plasmaFreq,
        omega,
        dep_cutoff,
        rayPowerAbsTol,
        dDeposition,
        maxSteps
    );

    cudaDeviceSynchronize();
}

// ************************************************************************* //
