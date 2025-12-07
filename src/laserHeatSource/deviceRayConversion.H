#ifndef DEVICE_RAY_CONVERSION_H
#define DEVICE_RAY_CONVERSION_H

#ifdef FOAM_USE_CUDA

#include "deviceRay.H"
#include "compactRay.H"

// Convert from OpenFOAM's compactRay to DeviceRay (host-side)
inline DeviceRay toDeviceRay(const Foam::compactRay& r)
{
    DeviceRay d;
    d.pos = make_double3(r.position_.x(), r.position_.y(), r.position_.z());
    d.dir = make_double3(r.direction_.x(), r.direction_.y(), r.direction_.z());
    d.power = r.power_;
    d.currentCell = r.currentCell_;
    d.globalIndex = r.globalRayIndex_;
    d.step = 0.0;
    return d;
}

// Convert from DeviceRay back to compactRay (host-side)
inline Foam::compactRay toCompactRay(const DeviceRay& d)
{
    return Foam::compactRay
    (
        Foam::point(d.pos.x, d.pos.y, d.pos.z),
        Foam::vector(d.dir.x, d.dir.y, d.dir.z),
        d.power
    );
}

#endif // FOAM_USE_CUDA

#endif // DEVICE_RAY_CONVERSION_H
