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

// This must match the layout of DeviceRay in laserHeatSourceGpu.C
struct DeviceRay
{
    double3 pos;
    double3 dir;
    double  power;
    int     currentCell;
    int     globalIndex;
    double  step;
};


// GPU ray kernal
__global__ void gpuRayStepKernel
(
    DeviceRay* rays,
    int nRays,
    const double* VI,
    int nCells,
    double3 bbMin,
    double3 bbMax,
    double rayPowerAbsTol
)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= nRays) return;

    DeviceRay& r = rays[i];

    // Cull by power first
    if (r.power <= rayPowerAbsTol)
    {
        r.power = 0.0;
        return;
    }

    const int c = r.currentCell;
    if (c < 0 || c >= nCells) return;

    const double V = VI[c];
    if (V <= 0.0) return;

    const double pi = 3.14159265358979323846;
    const double h  = cbrt(V);
    const double iterator_distance = (0.5/pi)*h;

    // Step along the ray direction
    r.pos.x += iterator_distance * r.dir.x;
    r.pos.y += iterator_distance * r.dir.y;
    r.pos.z += iterator_distance * r.dir.z;

    // Cull if we’ve left the global bounding box
    if
    (
        r.pos.x < bbMin.x || r.pos.x > bbMax.x ||
        r.pos.y < bbMin.y || r.pos.y > bbMax.y ||
        r.pos.z < bbMin.z || r.pos.z > bbMax.z
    )
    {
        r.power = 0.0;
    }
}


// C-linkage wrapper we can call from C++ code
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
)
{
    if (!dRays || nRays <= 0 || !dVI || nCells <= 0) return;

    const int blockSize = 128;
    const int gridSize  = (nRays + blockSize - 1)/blockSize;

    double3 bbMin = make_double3(bbMinX, bbMinY, bbMinZ);
    double3 bbMax = make_double3(bbMaxX, bbMaxY, bbMaxZ);

    gpuRayStepKernel<<<gridSize, blockSize>>>
    (
        dRays,
        nRays,
        dVI,
        nCells,
        bbMin,
        bbMax,
        rayPowerAbsTol
    );
    cudaDeviceSynchronize();
}


// ************************************************************************* //
