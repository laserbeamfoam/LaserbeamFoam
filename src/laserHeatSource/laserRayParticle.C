/*---------------------------------------------------------------------------*\
License
    This file is part of laserbeamFoam.

    laserbeamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    laserbeamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with laserbeamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "laserRayParticle.H"
#include "constants.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(laserRayParticle, 0);

    defineTemplateTypeNameAndDebugWithName
    (
        Cloud<laserRayParticle>,
        "laserRayParticleCloud",
        0
    );
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::laserRayParticle::laserRayParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label cellI,
    const vector& direction,
    const scalar power,
    const scalar dA,
    const label globalRayIndex
)
:
    particle(mesh, position, cellI),
    direction_(direction / (mag(direction) + VSMALL)),
    power_(power),
    initialPower_(power),
    dA_(dA),
    bounceCount_(0),
    globalRayIndex_(globalRayIndex),
    active_(true),
    path_()
{
    path_.append(position);
}


Foam::laserRayParticle::laserRayParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    particle(mesh, is, readFields, newFormat),
    direction_(vector::zero),
    power_(0),
    initialPower_(0),
    dA_(0),
    bounceCount_(0),
    globalRayIndex_(-1),
    active_(true),
    path_()
{
    if (readFields)
    {
        is  >> direction_
            >> power_
            >> initialPower_
            >> dA_
            >> bounceCount_
            >> globalRayIndex_
            >> active_;
    }

    is.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::laserRayParticle::move
(
    Cloud<laserRayParticle>& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    // Reference to the mesh
    const polyMesh& mesh = this->mesh();
    const scalarField& VI = mesh.cellVolumes();

    // References from tracking data
    const scalarField& alphaI = td.alphaFiltered_;
    const vectorField& nI = td.nFiltered_;

    // ================================================================
    // Main ray-tracing loop:
    //   - Physics check on current cell FIRST (interface/bulk)
    //   - Then track to the next face
    //   - Large displacement (mesh.bounds().mag()) so stepFraction
    //     advances slowly and the ray can cross many cells
    //   - After reflection, reset stepFraction to 0 so the new
    //     direction gets a fresh journey
    // ================================================================

    td.keepParticle = true;
    td.switchProcessor = false;

    // Ensure we have a start point. For newly created particles,
    // the constructor already appends the start position. For
    // particles arriving from another processor (deserialized via
    // Istream), path_ is empty — record the transfer position.
    if (path_.empty())
    {
        path_.append(position());
    }

    const scalar maxTrackLength = mesh.bounds().mag();

    while (td.keepParticle && !td.switchProcessor && active_)
    {
        const label cellI = cell();

        // Safety check
        if (cellI < 0 || cellI >= mesh.nCells())
        {
            active_ = false;
            td.keepParticle = false;
            break;
        }

        // Update visualisation fields
        td.rayQ_[cellI] += power_;
        td.rayNumber_[cellI] = globalRayIndex_;

        // Check if power is exhausted
        if (power_ < td.rayPowerAbsTol_ || power_ < SMALL)
        {
            path_.append(position());
            td.finishedRayIDs_.append(globalRayIndex_);
            td.finishedRayPaths_.append(path_);
            active_ = false;
            td.keepParticle = false;
            break;
        }

        // ---- Physics: check for interface or bulk metal ----

        if
        (
            mag(nI[cellI]) > 0.5
         && alphaI[cellI] >= td.depCutoff_
        )
        {
            // ================================================
            // Interface cell - Fresnel reflection/absorption
            // ================================================

            const scalar absorptivity = computeFresnelAbsorptivity
            (
                direction_,
                nI[cellI],
                td.plasmaFrequency_,
                td.angularFrequency_,
                td.resistivity_[cellI]
            );

            // Deposit absorbed fraction
            td.deposition_[cellI] += absorptivity * power_ / VI[cellI];

            // Reduce ray power
            power_ *= (1.0 - absorptivity);

            // Specular reflection about the interface normal
            vector n = nI[cellI];
            n /= mag(n);

            // Ensure normal points into the incident medium
            if ((direction_ & n) > 0)
            {
                n = -n;
            }

            direction_ = direction_ - 2.0*(direction_ & n)*n;
            direction_ /= mag(direction_) + VSMALL;

            bounceCount_++;
            path_.append(position());

            // If power is now below tolerance, kill the ray
            if (power_ < td.rayPowerAbsTol_)
            {
                td.finishedRayIDs_.append(globalRayIndex_);
                td.finishedRayPaths_.append(path_);
                active_ = false;
                td.keepParticle = false;
                break;
            }
        }
        else if (alphaI[cellI] >= td.depCutoff_)
        {
            // ================================================
            // Bulk metal - full absorption, ray dies
            // ================================================

            td.deposition_[cellI] += power_ / VI[cellI];
            power_ = 0.0;
            path_.append(position());
            td.finishedRayIDs_.append(globalRayIndex_);
            td.finishedRayPaths_.append(path_);
            active_ = false;
            td.keepParticle = false;
            break;
        }

        // ---- Track to the next face ----

        // Reset stepFraction before each tracking call. We are NOT
        // doing time-based tracking (rays are instantaneous), so
        // stepFraction should not accumulate across cell crossings.
        // Without this reset, stepFraction eventually reaches 1.0
        // and trackToAndHitFace stops advancing the particle, causing
        // an infinite loop that deadlocks parallel runs.
        stepFraction() = 0;
        const scalar f = 1.0;
        const vector displacement = direction_*maxTrackLength;

        trackToAndHitFace(displacement, f, cloud, td);
    }

    return td.keepParticle;
}


bool Foam::laserRayParticle::hitPatch
(
    Cloud<laserRayParticle>& cloud,
    trackingData& td
)
{
    // Return false: not handled here, proceed to specific patch handler
    // (hitProcessorPatch, hitWallPatch, etc.)
    return false;
}


void Foam::laserRayParticle::hitProcessorPatch
(
    Cloud<laserRayParticle>& cloud,
    trackingData& td
)
{
    // Store the path segment accumulated so far BEFORE the particle
    // is serialized and transferred (path_ is not serialized, so it
    // would be lost). Append current position as the end of this segment.
    path_.append(position());
    td.finishedRayIDs_.append(globalRayIndex_);
    td.finishedRayPaths_.append(path_);

    // The Cloud infrastructure will automatically transfer this
    // particle to the neighbouring processor. We just need to flag
    // the switch.
    td.switchProcessor = true;
}


void Foam::laserRayParticle::hitWallPatch
(
    Cloud<laserRayParticle>& cloud,
    trackingData& td
)
{
    // Wall boundary - absorb remaining power into the boundary cell
    // (or simply kill the ray if the wall is outside the material)
    const label cellI = cell();

    if (cellI >= 0 && power_ > SMALL)
    {
        const scalarField& VI = this->mesh().cellVolumes();
        td.deposition_[cellI] += power_ / VI[cellI];
        power_ = 0.0;
    }

    path_.append(position());
    td.finishedRayIDs_.append(globalRayIndex_);
    td.finishedRayPaths_.append(path_);
    active_ = false;
    td.keepParticle = false;
}


Foam::scalar Foam::laserRayParticle::computeFresnelAbsorptivity
(
    const vector& rayDirection,
    const vector& surfaceNormal,
    const scalar plasmaFrequency,
    const scalar angularFrequency,
    const scalar resistivity
)
{
    const scalar damping_frequency =
        plasmaFrequency*plasmaFrequency
       *constant::electromagnetic::epsilon0.value()
       *resistivity;

    const scalar e_r =
        1.0
      - (
            sqr(plasmaFrequency)
           /(sqr(angularFrequency) + sqr(damping_frequency))
        );

    const scalar e_i =
        (damping_frequency/angularFrequency)
       *(
            sqr(plasmaFrequency)
           /(sqr(angularFrequency) + sqr(damping_frequency))
        );

    const scalar ref_index =
        Foam::sqrt((Foam::sqrt(e_r*e_r + e_i*e_i) + e_r)/2.0);

    const scalar ext_coefficient =
        Foam::sqrt((Foam::sqrt(e_r*e_r + e_i*e_i) - e_r)/2.0);

    // Unit surface normal
    vector n = surfaceNormal;
    n /= mag(n);

    // Unit ray direction
    vector d = rayDirection;
    const scalar dMag = mag(d);
    if (dMag <= SMALL)
    {
        return 1.0;
    }
    d /= dMag;

    const vector kin = -d;
    scalar cosTheta = kin & n;

    if (cosTheta < 0.0)
    {
        n = -n;
        cosTheta = -cosTheta;
    }

    cosTheta = Foam::max(Foam::min(cosTheta, scalar(1.0)), scalar(0.0));
    const scalar theta_in = std::acos(cosTheta);
    const scalar sinTheta = Foam::sin(theta_in);
    const scalar cosTheta_in = Foam::cos(theta_in);

    const scalar alpha_laser = Foam::sqrt
    (
        Foam::sqrt
        (
            sqr(sqr(ref_index) - sqr(ext_coefficient) - sqr(sinTheta))
          + 4.0*sqr(ref_index)*sqr(ext_coefficient)
        )
      + sqr(ref_index)
      - sqr(ext_coefficient)
      - sqr(sinTheta)/2.0
    );

    const scalar beta_laser = Foam::sqrt
    (
        (
            Foam::sqrt
            (
                sqr(sqr(ref_index) - sqr(ext_coefficient) - sqr(sinTheta))
              + 4.0*sqr(ref_index)*sqr(ext_coefficient)
            )
          - sqr(ref_index)
          + sqr(ext_coefficient)
          + sqr(sinTheta)
        )/2.0
    );

    scalar R_s =
        (
            sqr(alpha_laser) + sqr(beta_laser)
          - 2.0*alpha_laser*cosTheta_in + sqr(cosTheta_in)
        )
       /(
            sqr(alpha_laser) + sqr(beta_laser)
          + 2.0*alpha_laser*cosTheta_in + sqr(cosTheta_in)
        );

    scalar R_p =
        R_s
       *(
            sqr(alpha_laser) + sqr(beta_laser)
          - 2.0*alpha_laser*sinTheta*Foam::tan(theta_in)
          + sqr(sinTheta)*sqr(Foam::tan(theta_in))
        )
       /(
            sqr(alpha_laser) + sqr(beta_laser)
          + 2.0*alpha_laser*sinTheta*Foam::tan(theta_in)
          + sqr(sinTheta)*sqr(Foam::tan(theta_in))
        );

    R_s = Foam::max(Foam::min(R_s, scalar(1.0)), scalar(0.0));
    R_p = Foam::max(Foam::min(R_p, scalar(1.0)), scalar(0.0));

    const scalar R = 0.5*(R_s + R_p);
    scalar absorptivity = 1.0 - R;

    return Foam::max(Foam::min(absorptivity, scalar(1.0)), scalar(0.0));
}


// * * * * * * * * * * * * * * * Read/Write  * * * * * * * * * * * * * * * * //

void Foam::laserRayParticle::readFields
(
    Cloud<laserRayParticle>& cloud
)
{
    bool valid = cloud.size();

    particle::readFields(cloud);

    IOField<vector> direction
    (
        cloud.newIOobject("direction", IOobject::MUST_READ),
        valid
    );
    cloud.checkFieldIOobject(cloud, direction);

    IOField<scalar> power
    (
        cloud.newIOobject("power", IOobject::MUST_READ),
        valid
    );
    cloud.checkFieldIOobject(cloud, power);

    IOField<scalar> initialPower
    (
        cloud.newIOobject("initialPower", IOobject::MUST_READ),
        valid
    );
    cloud.checkFieldIOobject(cloud, initialPower);

    IOField<scalar> dA
    (
        cloud.newIOobject("dA", IOobject::MUST_READ),
        valid
    );
    cloud.checkFieldIOobject(cloud, dA);

    IOField<label> bounceCount
    (
        cloud.newIOobject("bounceCount", IOobject::MUST_READ),
        valid
    );
    cloud.checkFieldIOobject(cloud, bounceCount);

    IOField<label> globalRayIndex
    (
        cloud.newIOobject("globalRayIndex", IOobject::MUST_READ),
        valid
    );
    cloud.checkFieldIOobject(cloud, globalRayIndex);

    label i = 0;
    for (laserRayParticle& p : cloud)
    {
        p.direction_ = direction[i];
        p.power_ = power[i];
        p.initialPower_ = initialPower[i];
        p.dA_ = dA[i];
        p.bounceCount_ = bounceCount[i];
        p.globalRayIndex_ = globalRayIndex[i];
        ++i;
    }
}


void Foam::laserRayParticle::writeFields
(
    const Cloud<laserRayParticle>& cloud
)
{
    particle::writeFields(cloud);

    const label np = cloud.size();

    IOField<vector> direction
    (
        cloud.newIOobject("direction", IOobject::NO_READ),
        np
    );
    IOField<scalar> power
    (
        cloud.newIOobject("power", IOobject::NO_READ),
        np
    );
    IOField<scalar> initialPower
    (
        cloud.newIOobject("initialPower", IOobject::NO_READ),
        np
    );
    IOField<scalar> dA
    (
        cloud.newIOobject("dA", IOobject::NO_READ),
        np
    );
    IOField<label> bounceCount
    (
        cloud.newIOobject("bounceCount", IOobject::NO_READ),
        np
    );
    IOField<label> globalRayIndex
    (
        cloud.newIOobject("globalRayIndex", IOobject::NO_READ),
        np
    );

    label i = 0;
    for (const laserRayParticle& p : cloud)
    {
        direction[i] = p.direction_;
        power[i] = p.power_;
        initialPower[i] = p.initialPower_;
        dA[i] = p.dA_;
        bounceCount[i] = p.bounceCount_;
        globalRayIndex[i] = p.globalRayIndex_;
        ++i;
    }

    direction.write(np > 0);
    power.write(np > 0);
    initialPower.write(np > 0);
    dA.write(np > 0);
    bounceCount.write(np > 0);
    globalRayIndex.write(np > 0);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const laserRayParticle& p
)
{
    os  << static_cast<const particle&>(p)
        << token::SPACE << p.direction_
        << token::SPACE << p.power_
        << token::SPACE << p.initialPower_
        << token::SPACE << p.dA_
        << token::SPACE << p.bounceCount_
        << token::SPACE << p.globalRayIndex_
        << token::SPACE << p.active_;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
