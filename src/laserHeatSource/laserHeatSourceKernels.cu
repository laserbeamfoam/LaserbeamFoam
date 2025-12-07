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

// laserHeatSourceKernels.cu
#include <cuda_runtime.h>

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

// A tiny no-op kernel: just proves we can launch something.
__global__ void noopRayKernel(DeviceRay* rays, int nRays)
{
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    if (i >= nRays) return;

    const double step = rays[i].step;

    rays[i].pos.x += step * rays[i].dir.x;
    rays[i].pos.y += step * rays[i].dir.y;
    rays[i].pos.z += step * rays[i].dir.z;
}

// C-linkage wrapper we can call from C++ code
extern "C"
void launchNoopRayKernel(DeviceRay* dRays, int nRays)
{
    if (!dRays || nRays <= 0) return;

    const int blockSize = 128;
    const int gridSize  = (nRays + blockSize - 1)/blockSize;

    noopRayKernel<<<gridSize, blockSize>>>(dRays, nRays);
    cudaDeviceSynchronize();
}

// ************************************************************************* //
