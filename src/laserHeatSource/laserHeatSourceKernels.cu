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

// =====================================================================
// device helper: voxelFindCellForRay (GPU equivalent of findCellVoxelised)
// =====================================================================
__device__ int voxelFindCellForRay
(
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
    const double*    resistivity,  // resistivity
    int              nCells,
    double3          bbMin,
    double3          bbMax,
    int3             voxelDims,
    double3          voxelSize,
    const int*       voxelCell,
    const double3*   cellMin,
    const double3*   cellMax,
    double           plasmaFreq,
    double           omega,
    double           dep_cutoff,
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
        const double3 nval  = n[cellI];

        // iterator distance like CPU
        const double step = (0.5/M_PI)*cbrt(VIc);

        // ------------------------------
        // Interface detection + reflection (no Fresnel yet)
        // ------------------------------
        const double nmag =
            sqrt(nval.x*nval.x + nval.y*nval.y + nval.z*nval.z);

        const bool atInterface = (nmag > 0.5 && alphac >= dep_cutoff);


        // =====================================================================
        // Fresnel interface absorption + reflection (GPU port of CPU logic)
        // =====================================================================
        if (atInterface)
        {
            // -------------------------------
            // 1. Unit normal (double3)
            // -------------------------------
            // Load the surface normal (double3) for this cell
            const double3 nc = n[cellI];
            const double invMagN =
                rsqrt(nc.x*nc.x + nc.y*nc.y + nc.z*nc.z + 1e-20);
            double Nx = nc.x * invMagN;
            double Ny = nc.y * invMagN;
            double Nz = nc.z * invMagN;

            // -------------------------------
            // 2. Normalise incoming direction
            // -------------------------------
            double Dx = r.dir.x;
            double Dy = r.dir.y;
            double Dz = r.dir.z;

            double invMagD = rsqrt(Dx*Dx + Dy*Dy + Dz*Dz + 1e-20);
            Dx *= invMagD;
            Dy *= invMagD;
            Dz *= invMagD;

            // Incoming direction is -d in CPU code
            double kinx = -Dx;
            double kiny = -Dy;
            double kinz = -Dz;

            // -------------------------------
            // 3. cos(theta) = kin · N
            // -------------------------------
            double cosTheta = kinx*Nx + kiny*Ny + kinz*Nz;

            // If negative, flip normal
            if (cosTheta < 0.0)
            {
                cosTheta = -cosTheta;
                // Flip normal
                // (we can just negate the components)
                // Nx,Ny,Nz are const after this block’s scope, so local copies:
                // but since Nx,Ny,Nz are locals, it's safe to modify
                Nx = -Nx;
                Ny = -Ny;
                Nz = -Nz;
            }

            // Clamp within [0,1]
            cosTheta = fmin(fmax(cosTheta, 0.0), 1.0);

            // theta_in
            const double theta_in = acos(cosTheta);
            const double sinTheta = sin(theta_in);

            // -------------------------------
            // 4. Material optical response
            // -------------------------------
            const double pf2 = plasmaFreq * plasmaFreq;
            const double om2 = omega * omega;

            const double epsF = 1e-12;
            const double damping =
                pf2 / (om2 + pf2) * (resistivity[cellI] * epsF);
            // CPU uses:
            // damping = plasmaFreq^2 * eps0 * resistivity

            const double eps_r =
                1.0 - (pf2 / (om2 + pf2));

            const double eps_i =
                (damping / omega) * (pf2 / (om2 + pf2));

            // refractive index real/imag parts
            const double root = sqrt(eps_r*eps_r + eps_i*eps_i);
            const double n_real = sqrt((root + eps_r) * 0.5);
            const double n_imag = sqrt((root - eps_r) * 0.5);

            // -------------------------------
            // 5. Compute Fresnel R_s and R_p
            // -------------------------------
            // CPU calls these alpha_laser and beta_laser:
            const double alpha_laser = n_real;
            const double beta_laser  = n_imag;

            // Using CPU formula:
            const double cosT = cosTheta;       // rename for clarity
            const double sinT = sinTheta;

            // R_s numerator/denominator
            const double Rs_num =
                alpha_laser*alpha_laser + beta_laser*beta_laser
                - 2.0*alpha_laser*cosT + cosT*cosT;

            const double Rs_den =
                alpha_laser*alpha_laser + beta_laser*beta_laser
                + 2.0*alpha_laser*cosT + cosT*cosT;

            double R_s = Rs_num / Rs_den;

            // R_p numerator/denominator
            const double Rp_num =
                R_s * (alpha_laser*alpha_laser + beta_laser*beta_laser
                       - 2.0*alpha_laser*sinT*cosT/sinT
                       + sinT*sinT*(cosT*cosT)/(sinT*sinT));  // simplified CPU expression

            const double Rp_den =
                alpha_laser*alpha_laser + beta_laser*beta_laser
                + 2.0*alpha_laser*sinT*cosT/sinT
                + sinT*sinT*(cosT*cosT)/(sinT*sinT);

            double R_p = Rp_num / Rp_den;

            // Clamp Fresnel reflectances
            R_s = fmin(fmax(R_s, 0.0), 1.0);
            R_p = fmin(fmax(R_p, 0.0), 1.0);

            const double R = 0.5*(R_s + R_p);
            double absorptivity = 1.0 - R;
            absorptivity = fmin(fmax(absorptivity, 0.0), 1.0);

            // -------------------------------
            // 6. Deposit absorbed energy
            // -------------------------------
            const double dQ = absorptivity * r.power;
            atomicAdd(&deposition[cellI], dQ / VIc);

            // Reduce ray power
            r.power -= dQ;
            if (r.power <= rayPowerAbsTol)
            {
                r.power = 0.0;
                break;
            }

            // -------------------------------
            // 7. Reflect direction
            // -------------------------------
            const double dDotN = Dx*Nx + Dy*Ny + Dz*Nz;

            double Rx = Dx - 2.0*dDotN*Nx;
            double Ry = Dy - 2.0*dDotN*Ny;
            double Rz = Dz - 2.0*dDotN*Nz;

            // Normalise reflected direction
            const double invMagR = rsqrt(Rx*Rx + Ry*Ry + Rz*Rz + 1e-20);
            r.dir.x = Rx * invMagR;
            r.dir.y = Ry * invMagR;
            r.dir.z = Rz * invMagR;

            // -------------------------------
            // 8. Push slightly forward
            // -------------------------------
            const double eps = 1e-6 * step;   // small nudge
            r.pos.x += eps * r.dir.x;
            r.pos.y += eps * r.dir.y;
            r.pos.z += eps * r.dir.z;

            continue;  // continue GPU ray marching loop
        }


        // ------------------------------
        // Bulk absorption (only if NOT interface)
        // ------------------------------
        if (!atInterface && alphac >= dep_cutoff)
        {
            const double volSrc = (VIc > 0.0 ? r.power/VIc : 0.0);

            atomicAdd(&deposition[cellI], volSrc);
            r.power = 0.0;
            break;  // ray is finished
        }

        // ------------------------------
        // Step the ray in its current direction (possibly reflected)
        // ------------------------------
        r.pos.x += step * r.dir.x;
        r.pos.y += step * r.dir.y;
        r.pos.z += step * r.dir.z;

        // count this iteration
        ++steps;

    }

    // after maxSteps or termination, the ray is done
}


// =====================================================================
// C-linkage launcher
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
