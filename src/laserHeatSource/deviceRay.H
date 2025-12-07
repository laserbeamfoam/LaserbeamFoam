#ifndef DEVICE_RAY_H
#define DEVICE_RAY_H

#include <cuda_runtime.h>

// Pure GPU-safe struct — no OpenFOAM headers
struct DeviceRay
{
    double3 pos;
    double3 dir;
    double  power;
    int     currentCell;
    int     globalIndex;
    double  step;

    __host__ __device__
    DeviceRay() {}
};

#endif // DEVICE_RAY_H
