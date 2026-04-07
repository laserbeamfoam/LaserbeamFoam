/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 LaserbeamFoam contributors
-------------------------------------------------------------------------------
License
    This file is part of laserbeamFoam.

    laserbeamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    laserbeamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with laserbeamFoam. If not, see <http://www.gnu.org/licenses/>.

Description
    Implementation of the laser heat-source profile classes.

\*---------------------------------------------------------------------------*/

#include "HeatSourceProfiles.H"

#include "IFstream.H"
#include "StringStream.H"
#include "stringOps.H"
#include "unitConversion.H"
#include "wordList.H"

namespace Foam
{

namespace
{
    SamplingStrategy readSamplingStrategy(const word& rawType)
    {
        const word type(word(stringOps::lower(rawType)));

        if (type == "radialpolar" || type == "radial_polar")
        {
            return SamplingStrategy::RADIAL_POLAR;
        }
        else if (type == "cartesian")
        {
            return SamplingStrategy::CARTESIAN;
        }

        FatalErrorInFunction
            << "Unknown profile sampling type '" << rawType << "'. "
            << "Expected 'radialPolar' or 'cartesian'."
            << exit(FatalError);

        return SamplingStrategy::RADIAL_POLAR;
    }


    void beamLocalBasis
    (
        const dictionary& laserDict,
        vector& u,
        vector& v,
        vector& beamAxis
    )
    {
        beamAxis =
            laserDict.lookupOrDefault<vector>("V_incident", vector(0, 1, 0));
        beamAxis /= (mag(beamAxis) + VSMALL);

        const vector a =
            (mag(beamAxis.z()) < 0.9) ? vector(0, 0, 1) : vector(0, 1, 0);

        u = beamAxis ^ a;
        u /= (mag(u) + VSMALL);
        v = beamAxis ^ u;
        v /= (mag(v) + VSMALL);
    }


    scalar defaultProfileScale
    (
        const dictionary& profileDict,
        const dictionary& laserDict,
        const scalar fallback
    )
    {
        scalar scale = profileDict.lookupOrDefault<scalar>("radius", scalar(-1));

        if (scale <= SMALL)
        {
            scale = profileDict.lookupOrDefault<scalar>
            (
                "laserRadius",
                laserDict.lookupOrDefault<scalar>("laserRadius", scalar(-1))
            );
        }

        if (scale <= SMALL)
        {
            scale = fallback;
        }

        return max(scale, SMALL);
    }


    vector defaultCartesianHalfWidths
    (
        const dictionary& dict,
        const vector& fallback
    )
    {
        scalar halfWidth =
            dict.lookupOrDefault<scalar>
            (
                "halfWidth",
                dict.lookupOrDefault<scalar>
                (
                    "xHalfWidth",
                    dict.lookupOrDefault<scalar>
                    (
                        "xLength",
                        scalar(2)*fallback.x()
                    )/scalar(2)
                )
            );

        scalar halfHeight =
            dict.lookupOrDefault<scalar>
            (
                "halfHeight",
                dict.lookupOrDefault<scalar>
                (
                    "zHalfWidth",
                    dict.lookupOrDefault<scalar>
                    (
                        "zLength",
                        scalar(2)*fallback.z()
                    )/scalar(2)
                )
            );

        halfWidth = max(halfWidth, SMALL);
        halfHeight = max(halfHeight, SMALL);

        return vector(halfWidth, 0, halfHeight);
    }


    scalar blendFactor(const scalar x)
    {
        const scalar clamped = min(max(x, scalar(0)), scalar(1));
        return clamped*clamped*(scalar(3) - scalar(2)*clamped);
    }


    void readSymmetric2x2
    (
        const dictionary& dict,
        scalar& xx,
        scalar& xz,
        scalar& zz
    )
    {
        xx = dict.get<scalar>("xx");
        xz = dict.lookupOrDefault<scalar>("xz", scalar(0));
        zz = dict.get<scalar>("zz");
    }


    scalar smallestEigenvalue2x2
    (
        const scalar xx,
        const scalar xz,
        const scalar zz
    )
    {
        const scalar trace = xx + zz;
        const scalar disc =
            Foam::sqrt(max(sqr(xx - zz) + scalar(4)*sqr(xz), scalar(0)));

        return scalar(0.5)*(trace - disc);
    }


    void invertSymmetric2x2
    (
        const scalar xx,
        const scalar xz,
        const scalar zz,
        scalar& invXX,
        scalar& invXZ,
        scalar& invZZ
    )
    {
        const scalar det = xx*zz - sqr(xz);

        if (mag(det) <= SMALL)
        {
            FatalErrorInFunction
                << "The supplied 2x2 matrix is singular and cannot be inverted."
                << exit(FatalError);
        }

        invXX = zz/det;
        invXZ = -xz/det;
        invZZ = xx/det;
    }


    scalar besselFunction
    (
        const label order,
        const scalar value
    )
    {
        if (order == 0)
        {
            return j0(value);
        }
        else if (order == 1)
        {
            return j1(value);
        }

        return jn(order, value);
    }


    scalar estimateBesselZero
    (
        const label azimuthalOrder,
        const label radialOrder
    )
    {
        static const scalar known[4][4] =
        {
            {2.4048255577, 5.5200781103, 8.6537279129, 11.7915344391},
            {3.8317059702, 7.0155866698, 10.1734681351, 13.3236919363},
            {5.1356223018, 8.4172441404, 11.6198411721, 14.7959517824},
            {6.3801618952, 9.7610231299, 13.0152007217, 16.2234661603}
        };

        if
        (
            azimuthalOrder >= 0 && azimuthalOrder < 4
         && radialOrder >= 1 && radialOrder <= 4
        )
        {
            return known[azimuthalOrder][radialOrder - 1];
        }

        return constant::mathematical::pi
            *(scalar(radialOrder) + scalar(0.5)*scalar(azimuthalOrder) - 0.25);
    }


    autoPtr<HeatSourceProfileModifier> makeHeatSourceProfileModifier
    (
        const word& rawType,
        const dictionary& dict,
        const dictionary& laserDict,
        const scalar defaultScale,
        const vector& cartesianHalfWidths
    )
    {
        word type(rawType);
        type = word(stringOps::lower(type));

        if (type == "planebias" || type == "tiltbias")
        {
            return autoPtr<HeatSourceProfileModifier>
            (
                new PlaneBiasModifier(dict, laserDict, defaultScale)
            );
        }
        else if (type == "drumnoise" || type == "besselnoise")
        {
            return autoPtr<HeatSourceProfileModifier>
            (
                new DrumNoiseModifier(dict, defaultScale)
            );
        }
        else if (type == "cartesianmodenoise")
        {
            return autoPtr<HeatSourceProfileModifier>
            (
                new CartesianModeNoiseModifier(dict, cartesianHalfWidths)
            );
        }
        else if (type == "gridnoise")
        {
            return autoPtr<HeatSourceProfileModifier>
            (
                new GridNoiseModifier(dict, cartesianHalfWidths)
            );
        }

        FatalErrorInFunction
            << "Unknown heat source profile modifier type '" << rawType << "'"
            << nl << dict << exit(FatalError);

        return autoPtr<HeatSourceProfileModifier>();
    }
}


DrumNoiseMode::DrumNoiseMode()
:
    azimuthalOrder(0),
    radialOrder(1),
    amplitude(0),
    phase(0),
    frequency(0),
    currentAmplitude(0),
    currentPhase(0),
    normalisation(1),
    meanValue(0)
{}


CartesianNoiseMode::CartesianNoiseMode()
:
    xOrder(1),
    zOrder(0),
    amplitude(0),
    phase(0),
    frequency(0),
    currentAmplitude(0),
    currentPhase(0)
{}


HeatSourceProfileSampling::HeatSourceProfileSampling()
:
    type_(SamplingStrategy::RADIAL_POLAR),
    radius_(0),
    nRadial_(0),
    nAngular_(0),
    halfWidthX_(0),
    halfWidthZ_(0),
    nCartesianX_(0),
    nCartesianZ_(0),
    minRayPowerFraction_(0)
{}


HeatSourceProfileSampling::HeatSourceProfileSampling
(
    const dictionary& profileDict
)
:
    HeatSourceProfileSampling()
{
    if (!profileDict.found("sampling"))
    {
        FatalErrorInFunction
            << "Each heat source profile must define a 'sampling' sub-dictionary."
            << nl << profileDict << exit(FatalError);
    }

    const dictionary& samplingDict = profileDict.subDict("sampling");
    const word rawType = samplingDict.get<word>("type");
    type_ = readSamplingStrategy(rawType);
    minRayPowerFraction_ =
        samplingDict.lookupOrDefault<scalar>("minRayPowerFraction", scalar(0));

    if (minRayPowerFraction_ < 0 || minRayPowerFraction_ >= 1)
    {
        FatalErrorInFunction
            << "sampling.minRayPowerFraction must satisfy 0 <= value < 1."
            << nl << samplingDict << exit(FatalError);
    }

    if (type_ == SamplingStrategy::RADIAL_POLAR)
    {
        radius_ = samplingDict.get<scalar>("radius");
        nRadial_ = samplingDict.get<label>("nRadial");
        nAngular_ = samplingDict.get<label>("nAngular");

        if (radius_ <= SMALL || nRadial_ < 1 || nAngular_ < 1)
        {
            FatalErrorInFunction
                << "Radial-polar sampling requires radius > 0 and "
                << "positive nRadial/nAngular values."
                << nl << samplingDict << exit(FatalError);
        }
    }
    else
    {
        halfWidthX_ = samplingDict.get<scalar>("halfWidthX");
        halfWidthZ_ = samplingDict.get<scalar>("halfWidthZ");
        nCartesianX_ = samplingDict.get<label>("nCartesianX");
        nCartesianZ_ = samplingDict.get<label>("nCartesianZ");

        if
        (
            halfWidthX_ <= SMALL
         || halfWidthZ_ <= SMALL
         || nCartesianX_ < 1
         || nCartesianZ_ < 1
        )
        {
            FatalErrorInFunction
                << "Cartesian sampling requires positive half widths and "
                << "positive nCartesianX/nCartesianZ values."
                << nl << samplingDict << exit(FatalError);
        }
    }
}


PlaneBiasModifier::PlaneBiasModifier
(
    const dictionary& dict,
    const dictionary& laserDict,
    const scalar defaultScale
)
:
    localDirection_(vector::zero),
    localOrigin_
    (
        dict.lookupOrDefault<vector>
        (
            "originLocal",
            dict.lookupOrDefault<vector>("centerShift", vector::zero)
        )
    ),
    amplitude_(dict.lookupOrDefault<scalar>("amplitude", scalar(0))),
    offset_(dict.lookupOrDefault<scalar>("offset", scalar(0))),
    referenceScale_
    (
        max
        (
            dict.lookupOrDefault<scalar>("referenceScale", defaultScale),
            SMALL
        )
    ),
    minFactor_(dict.lookupOrDefault<scalar>("minFactor", scalar(0)))
{
    vector u(vector::zero);
    vector v(vector::zero);
    vector beamAxis(vector::zero);
    beamLocalBasis(laserDict, u, v, beamAxis);

    vector planeNormal =
        dict.lookupOrDefault<vector>("planeNormal", beamAxis);
    planeNormal /= (mag(planeNormal) + VSMALL);

    vector tiltDirection =
        dict.lookupOrDefault<vector>
        (
            "tiltDirection",
            dict.lookupOrDefault<vector>("direction", u)
        );

    tiltDirection -= (tiltDirection & planeNormal)*planeNormal;
    tiltDirection -= (tiltDirection & beamAxis)*beamAxis;

    if (mag(tiltDirection) <= SMALL)
    {
        tiltDirection = u;
    }

    tiltDirection /= (mag(tiltDirection) + VSMALL);

    localDirection_ =
        vector(tiltDirection & u, scalar(0), tiltDirection & v);

    if (mag(localDirection_) <= SMALL)
    {
        localDirection_ = vector(1, 0, 0);
    }
    else
    {
        localDirection_ /= (mag(localDirection_) + VSMALL);
    }
}


scalar PlaneBiasModifier::evaluate
(
    const vector& localPosition,
    const scalar
) const
{
    const vector shifted = localPosition - localOrigin_;
    const scalar coord =
        (shifted.x()*localDirection_.x() + shifted.z()*localDirection_.z())
       /referenceScale_;

    return max(scalar(1) + offset_ + amplitude_*coord, minFactor_);
}


DrumNoiseModifier::DrumNoiseModifier
(
    const dictionary& dict,
    const scalar defaultRadius
)
:
    modes_(0),
    radius_(max(dict.lookupOrDefault<scalar>("radius", defaultRadius), SMALL)),
    minFactor_(dict.lookupOrDefault<scalar>("minFactor", scalar(0))),
    randomPhase_
    (
        dict.lookupOrDefault<Switch>
        (
            "randomPhase",
            dict.lookupOrDefault<bool>("stepwise", false)
        )
    ),
    randomAmplitude_(dict.lookupOrDefault<Switch>("randomAmplitude", false)),
    amplitudeJitter_(max(dict.lookupOrDefault<scalar>("amplitudeJitter", scalar(0)), scalar(0))),
    evolution_
    (
        word
        (
            stringOps::lower
            (
                dict.lookupOrDefault<word>
                (
                    "evolution",
                    dict.lookupOrDefault<bool>("stepwise", false)
                  ? word("stepwise")
                  : word("continuous")
                )
            )
        )
    ),
    lastUpdateTime_(-GREAT),
    rndGen_(dict.lookupOrDefault<label>("seed", Random::defaultSeed))
{
    if (evolution_ != "continuous" && evolution_ != "stepwise")
    {
        FatalErrorInFunction
            << "DrumNoiseModifier expects evolution to be either "
            << "'continuous' or 'stepwise'."
            << nl << dict << exit(FatalError);
    }

    if (!dict.found("modes"))
    {
        FatalErrorInFunction
            << "DrumNoiseModifier requires a non-empty 'modes' dictionary."
            << nl << dict << exit(FatalError);
    }

    const dictionary& modesDict = dict.subDict("modes");
    const wordList keys = modesDict.toc();

    if (keys.empty())
    {
        FatalErrorInFunction
            << "DrumNoiseModifier requires at least one mode."
            << nl << dict << exit(FatalError);
    }

    modes_.setSize(keys.size());

    forAll(keys, i)
    {
        const dictionary& modeDict = modesDict.subDict(keys[i]);
        DrumNoiseMode& mode = modes_[i];

        mode.azimuthalOrder =
            modeDict.lookupOrDefault<label>("azimuthalOrder", 0);
        mode.radialOrder =
            modeDict.lookupOrDefault<label>("radialOrder", 1);
        mode.amplitude =
            modeDict.lookupOrDefault<scalar>
            (
                "amplitude",
                modeDict.lookupOrDefault<scalar>("weight", scalar(0))
            );
        mode.phase = degToRad(modeDict.lookupOrDefault<scalar>("phase", scalar(0)));
        mode.frequency =
            modeDict.lookupOrDefault<scalar>("frequency", scalar(0));
        mode.currentAmplitude = mode.amplitude;
        mode.currentPhase = mode.phase;

        if (mode.radialOrder < 1 || mode.azimuthalOrder < 0)
        {
            FatalErrorInFunction
                << "Drum noise modes require azimuthalOrder >= 0 and "
                << "radialOrder >= 1."
                << nl << modeDict << exit(FatalError);
        }

        const scalar alpha =
            estimateBesselZero(mode.azimuthalOrder, mode.radialOrder);

        scalar maxAbsValue = SMALL;
        scalar meanValue = 0.0;
        const label nSamples = 512;

        for (label j = 0; j <= nSamples; ++j)
        {
            const scalar xi = scalar(j)/scalar(nSamples);
            const scalar value =
                besselFunction(mode.azimuthalOrder, alpha*xi);

            maxAbsValue = max(maxAbsValue, mag(value));

            if (mode.azimuthalOrder == 0)
            {
                const scalar w =
                    (j == 0 || j == nSamples) ? scalar(0.5) : scalar(1);
                meanValue += w*value*xi;
            }
        }

        mode.normalisation = max(maxAbsValue, SMALL);
        mode.meanValue =
            mode.azimuthalOrder == 0
          ? scalar(2)*meanValue/scalar(nSamples)/mode.normalisation
          : scalar(0);
    }
}


scalar DrumNoiseModifier::modeSignal
(
    const DrumNoiseMode& mode,
    const scalar r,
    const scalar theta
) const
{
    const scalar alpha =
        estimateBesselZero(mode.azimuthalOrder, mode.radialOrder);
    const scalar radialValue =
        besselFunction(mode.azimuthalOrder, alpha*r/radius_)
       /mode.normalisation;
    const scalar angularValue =
        mode.azimuthalOrder == 0
      ? scalar(1)
      : Foam::cos
        (
            scalar(mode.azimuthalOrder)*theta + mode.currentPhase
        );

    scalar value = radialValue*angularValue;

    if (mode.azimuthalOrder == 0)
    {
        value -= mode.meanValue;
    }

    return mode.currentAmplitude*value;
}


scalar DrumNoiseModifier::evaluate
(
    const vector& localPosition,
    const scalar
) const
{
    const scalar r =
        Foam::sqrt(sqr(localPosition.x()) + sqr(localPosition.z()));

    if (r > radius_)
    {
        return scalar(1);
    }

    const scalar theta =
        (r > SMALL)
      ? Foam::atan2(localPosition.z(), localPosition.x())
      : scalar(0);

    scalar modulation = 0.0;

    forAll(modes_, i)
    {
        modulation += modeSignal(modes_[i], r, theta);
    }

    return max(scalar(1) + modulation, minFactor_);
}


void DrumNoiseModifier::update(const scalar currentTime, const scalar)
{
    if (evolution_ == "stepwise")
    {
        if (mag(currentTime - lastUpdateTime_) < SMALL)
        {
            return;
        }

        forAll(modes_, i)
        {
            DrumNoiseMode& mode = modes_[i];

            mode.currentPhase =
                randomPhase_
              ? constant::mathematical::twoPi*rndGen_.sample01<scalar>()
              : mode.phase;

            if (randomAmplitude_)
            {
                const scalar jitter =
                    amplitudeJitter_*rndGen_.GaussNormal<scalar>();
                mode.currentAmplitude = mode.amplitude*(scalar(1) + jitter);
            }
            else
            {
                mode.currentAmplitude = mode.amplitude;
            }
        }
    }
    else
    {
        forAll(modes_, i)
        {
            DrumNoiseMode& mode = modes_[i];
            const scalar temporalPhase =
                constant::mathematical::twoPi*mode.frequency*currentTime;

            mode.currentPhase = mode.phase + temporalPhase;
            mode.currentAmplitude =
                mode.amplitude*Foam::cos(temporalPhase);
        }
    }

    lastUpdateTime_ = currentTime;
}


CartesianModeNoiseModifier::CartesianModeNoiseModifier
(
    const dictionary& dict,
    const vector& defaultHalfWidths
)
:
    modes_(0),
    halfWidth_(defaultCartesianHalfWidths(dict, defaultHalfWidths).x()),
    halfHeight_(defaultCartesianHalfWidths(dict, defaultHalfWidths).z()),
    minFactor_(dict.lookupOrDefault<scalar>("minFactor", scalar(0))),
    randomPhase_
    (
        dict.lookupOrDefault<Switch>
        (
            "randomPhase",
            dict.lookupOrDefault<bool>("stepwise", false)
        )
    ),
    randomAmplitude_(dict.lookupOrDefault<Switch>("randomAmplitude", false)),
    amplitudeJitter_
    (
        max(dict.lookupOrDefault<scalar>("amplitudeJitter", scalar(0)), scalar(0))
    ),
    evolution_
    (
        word
        (
            stringOps::lower
            (
                dict.lookupOrDefault<word>
                (
                    "evolution",
                    dict.lookupOrDefault<bool>("stepwise", false)
                  ? word("stepwise")
                  : word("continuous")
                )
            )
        )
    ),
    lastUpdateTime_(-GREAT),
    rndGen_(dict.lookupOrDefault<label>("seed", Random::defaultSeed))
{
    if (evolution_ != "continuous" && evolution_ != "stepwise")
    {
        FatalErrorInFunction
            << "CartesianModeNoiseModifier expects evolution to be either "
            << "'continuous' or 'stepwise'."
            << nl << dict << exit(FatalError);
    }

    if (!dict.found("modes"))
    {
        FatalErrorInFunction
            << "CartesianModeNoiseModifier requires a non-empty 'modes' "
            << "dictionary."
            << nl << dict << exit(FatalError);
    }

    const dictionary& modesDict = dict.subDict("modes");
    const wordList keys = modesDict.toc();

    if (keys.empty())
    {
        FatalErrorInFunction
            << "CartesianModeNoiseModifier requires at least one mode."
            << nl << dict << exit(FatalError);
    }

    modes_.setSize(keys.size());

    forAll(keys, i)
    {
        const dictionary& modeDict = modesDict.subDict(keys[i]);
        CartesianNoiseMode& mode = modes_[i];

        mode.xOrder = modeDict.lookupOrDefault<label>("xOrder", 0);
        mode.zOrder = modeDict.lookupOrDefault<label>("zOrder", 0);
        mode.amplitude =
            modeDict.lookupOrDefault<scalar>
            (
                "amplitude",
                modeDict.lookupOrDefault<scalar>("weight", scalar(0))
            );
        mode.phase =
            degToRad(modeDict.lookupOrDefault<scalar>("phase", scalar(0)));
        mode.frequency =
            modeDict.lookupOrDefault<scalar>("frequency", scalar(0));
        mode.currentAmplitude = mode.amplitude;
        mode.currentPhase = mode.phase;

        if (mode.xOrder < 0 || mode.zOrder < 0)
        {
            FatalErrorInFunction
                << "Cartesian noise modes require xOrder >= 0 and zOrder >= 0."
                << nl << modeDict << exit(FatalError);
        }

        if (mode.xOrder == 0 && mode.zOrder == 0)
        {
            FatalErrorInFunction
                << "Cartesian noise modes cannot use xOrder = 0 and zOrder = 0 "
                << "because that would create a constant bias instead of "
                << "redistributing intensity."
                << nl << modeDict << exit(FatalError);
        }
    }
}


scalar CartesianModeNoiseModifier::modeSignal
(
    const CartesianNoiseMode& mode,
    const scalar x,
    const scalar z
) const
{
    const scalar xi =
        min(max((x + halfWidth_)/(scalar(2)*halfWidth_), scalar(0)), scalar(1));
    const scalar zi =
        min(max((z + halfHeight_)/(scalar(2)*halfHeight_), scalar(0)), scalar(1));

    const scalar xBasis =
        mode.xOrder == 0
      ? scalar(1)
      : Foam::cos
        (
            constant::mathematical::twoPi*scalar(mode.xOrder)*xi
          + mode.currentPhase
        );

    const scalar zBasis =
        mode.zOrder == 0
      ? scalar(1)
      : Foam::cos
        (
            constant::mathematical::twoPi*scalar(mode.zOrder)*zi
          + mode.currentPhase
        );

    return mode.currentAmplitude*xBasis*zBasis;
}


scalar CartesianModeNoiseModifier::evaluate
(
    const vector& localPosition,
    const scalar
) const
{
    if
    (
        mag(localPosition.x()) > halfWidth_ + SMALL
     || mag(localPosition.z()) > halfHeight_ + SMALL
    )
    {
        return scalar(1);
    }

    scalar modulation = 0.0;

    forAll(modes_, i)
    {
        modulation += modeSignal(modes_[i], localPosition.x(), localPosition.z());
    }

    return max(scalar(1) + modulation, minFactor_);
}


void CartesianModeNoiseModifier::update
(
    const scalar currentTime,
    const scalar
)
{
    if (evolution_ == "stepwise")
    {
        if (mag(currentTime - lastUpdateTime_) < SMALL)
        {
            return;
        }

        forAll(modes_, i)
        {
            CartesianNoiseMode& mode = modes_[i];

            mode.currentPhase =
                randomPhase_
              ? constant::mathematical::twoPi*rndGen_.sample01<scalar>()
              : mode.phase;

            if (randomAmplitude_)
            {
                const scalar jitter =
                    amplitudeJitter_*rndGen_.GaussNormal<scalar>();
                mode.currentAmplitude = mode.amplitude*(scalar(1) + jitter);
            }
            else
            {
                mode.currentAmplitude = mode.amplitude;
            }
        }
    }
    else
    {
        forAll(modes_, i)
        {
            CartesianNoiseMode& mode = modes_[i];

            if (lastUpdateTime_ < -scalar(0.5)*GREAT)
            {
                if (randomPhase_)
                {
                    mode.phase =
                        constant::mathematical::twoPi*rndGen_.sample01<scalar>();
                }

                if (randomAmplitude_)
                {
                    const scalar jitter =
                        amplitudeJitter_*rndGen_.GaussNormal<scalar>();
                    mode.amplitude *= scalar(1) + jitter;
                }
            }

            const scalar temporalPhase =
                constant::mathematical::twoPi*mode.frequency*currentTime;

            mode.currentPhase = mode.phase + temporalPhase;
            mode.currentAmplitude = mode.amplitude*Foam::cos(temporalPhase);
        }
    }

    lastUpdateTime_ = currentTime;
}


GridNoiseModifier::GridNoiseModifier
(
    const dictionary& dict,
    const vector& defaultHalfWidths
)
:
    nx_(dict.lookupOrDefault<label>("nx", 33)),
    nz_(dict.lookupOrDefault<label>("nz", dict.lookupOrDefault<label>("ny", 33))),
    halfWidth_(defaultCartesianHalfWidths(dict, defaultHalfWidths).x()),
    halfHeight_(defaultCartesianHalfWidths(dict, defaultHalfWidths).z()),
    amplitude_(max(dict.lookupOrDefault<scalar>("amplitude", scalar(0)), scalar(0))),
    minFactor_(dict.lookupOrDefault<scalar>("minFactor", scalar(0))),
    correlationLength_
    (
        max(dict.lookupOrDefault<scalar>("correlationLength", scalar(0)), scalar(0))
    ),
    smoothPasses_(max(dict.lookupOrDefault<label>("smoothPasses", 0), 0)),
    frequency_(max(dict.lookupOrDefault<scalar>("frequency", scalar(0)), scalar(0))),
    evolution_
    (
        word
        (
            stringOps::lower
            (
                dict.lookupOrDefault<word>
                (
                    "evolution",
                    dict.lookupOrDefault<bool>("stepwise", false)
                  ? word("stepwise")
                  : word("continuous")
                )
            )
        )
    ),
    currentBlend_(0),
    anchorTime_(-GREAT),
    nextAnchorTime_(-GREAT),
    lastUpdateTime_(-GREAT),
    currentMask_(0),
    nextMask_(0),
    rndGen_(dict.lookupOrDefault<label>("seed", Random::defaultSeed))
{
    if (nx_ <= 1 || nz_ <= 1)
    {
        FatalErrorInFunction
            << "GridNoiseModifier requires nx > 1 and nz > 1."
            << nl << dict << exit(FatalError);
    }

    if (evolution_ != "continuous" && evolution_ != "stepwise")
    {
        FatalErrorInFunction
            << "GridNoiseModifier expects evolution to be either "
            << "'continuous' or 'stepwise'."
            << nl << dict << exit(FatalError);
    }

    currentMask_.setSize(nx_*nz_, scalar(0));
    nextMask_.setSize(nx_*nz_, scalar(0));
}


void GridNoiseModifier::generateMask(scalarField& mask)
{
    mask.setSize(nx_*nz_, scalar(0));

    scalar sum = 0.0;

    forAll(mask, i)
    {
        const scalar value = rndGen_.GaussNormal<scalar>();
        mask[i] = value;
        sum += value;
    }

    const scalar mean = sum/scalar(mask.size());

    forAll(mask, i)
    {
        mask[i] -= mean;
    }

    label passes = smoothPasses_;

    if (correlationLength_ > SMALL)
    {
        const scalar dx = scalar(2)*halfWidth_/scalar(max(nx_ - 1, 1));
        const scalar dz = scalar(2)*halfHeight_/scalar(max(nz_ - 1, 1));
        const scalar h = max(dx, dz);
        passes =
            max
            (
                passes,
                label(correlationLength_/max(h, SMALL) + scalar(0.5))
            );
    }

    for (label pass = 0; pass < passes; ++pass)
    {
        scalarField blurred(mask.size(), scalar(0));

        for (label iz = 0; iz < nz_; ++iz)
        {
            for (label ix = 0; ix < nx_; ++ix)
            {
                scalar localSum = 0.0;
                label count = 0;

                for (label j = max(iz - 1, 0); j <= min(iz + 1, nz_ - 1); ++j)
                {
                    for (label i = max(ix - 1, 0); i <= min(ix + 1, nx_ - 1); ++i)
                    {
                        localSum += mask[j*nx_ + i];
                        ++count;
                    }
                }

                blurred[iz*nx_ + ix] = localSum/scalar(count);
            }
        }

        mask.transfer(blurred);
    }

    sum = 0.0;
    scalar maxAbs = SMALL;

    forAll(mask, i)
    {
        sum += mask[i];
    }

    const scalar blurredMean = sum/scalar(mask.size());

    forAll(mask, i)
    {
        mask[i] -= blurredMean;
        maxAbs = max(maxAbs, mag(mask[i]));
    }

    forAll(mask, i)
    {
        mask[i] /= maxAbs;
    }
}


scalar GridNoiseModifier::sampleMask
(
    const scalarField& mask,
    const scalar x,
    const scalar z
) const
{
    scalar fx = (x + halfWidth_)/(scalar(2)*halfWidth_)*scalar(nx_ - 1);
    scalar fz = (z + halfHeight_)/(scalar(2)*halfHeight_)*scalar(nz_ - 1);

    fx = min(max(fx, scalar(0)), scalar(nx_ - 1));
    fz = min(max(fz, scalar(0)), scalar(nz_ - 1));

    const label ix0 = min(label(fx), nx_ - 1);
    const label iz0 = min(label(fz), nz_ - 1);
    const label ix1 = min(ix0 + 1, nx_ - 1);
    const label iz1 = min(iz0 + 1, nz_ - 1);
    const scalar tx = fx - scalar(ix0);
    const scalar tz = fz - scalar(iz0);

    const scalar v00 = mask[iz0*nx_ + ix0];
    const scalar v10 = mask[iz0*nx_ + ix1];
    const scalar v01 = mask[iz1*nx_ + ix0];
    const scalar v11 = mask[iz1*nx_ + ix1];
    const scalar interpX0 = v00 + tx*(v10 - v00);
    const scalar interpX1 = v01 + tx*(v11 - v01);

    return interpX0 + tz*(interpX1 - interpX0);
}


scalar GridNoiseModifier::evaluate
(
    const vector& localPosition,
    const scalar
) const
{
    if
    (
        mag(localPosition.x()) > halfWidth_ + SMALL
     || mag(localPosition.z()) > halfHeight_ + SMALL
    )
    {
        return scalar(1);
    }

    scalar maskValue = sampleMask(currentMask_, localPosition.x(), localPosition.z());

    if (evolution_ == "continuous")
    {
        const scalar nextValue =
            sampleMask(nextMask_, localPosition.x(), localPosition.z());
        maskValue += currentBlend_*(nextValue - maskValue);
    }

    return max(scalar(1) + amplitude_*maskValue, minFactor_);
}


void GridNoiseModifier::update(const scalar currentTime, const scalar deltaT)
{
    if (evolution_ == "stepwise")
    {
        if (mag(currentTime - lastUpdateTime_) < SMALL)
        {
            return;
        }

        generateMask(currentMask_);
        lastUpdateTime_ = currentTime;
        return;
    }

    scalar interval =
        frequency_ > SMALL
      ? scalar(1)/frequency_
      : max(deltaT, SMALL);

    if (anchorTime_ < -scalar(0.5)*GREAT)
    {
        generateMask(currentMask_);
        generateMask(nextMask_);
        anchorTime_ = currentTime;
        nextAnchorTime_ = currentTime + interval;
    }

    while (currentTime >= nextAnchorTime_ - SMALL)
    {
        currentMask_ = nextMask_;
        generateMask(nextMask_);
        anchorTime_ = nextAnchorTime_;
        nextAnchorTime_ += interval;
    }

    currentBlend_ =
        blendFactor
        (
            (currentTime - anchorTime_)
           /max(nextAnchorTime_ - anchorTime_, SMALL)
        );
    lastUpdateTime_ = currentTime;
}


HeatSourceProfile::HeatSourceProfile(const dictionary& dict)
:
    weight_(dict.lookupOrDefault<scalar>("weight", 1.0)),
    sampling_(dict),
    modifiers_(0)
{
    if (weight_ < 0)
    {
        FatalErrorInFunction
            << "Heat source profile weights must be non-negative."
            << nl << dict << exit(FatalError);
    }
}


void HeatSourceProfile::buildModifiers
(
    const dictionary& profileDict,
    const dictionary& laserDict,
    const scalar defaultScale,
    const vector& cartesianHalfWidths
)
{
    if (!profileDict.found("modifiers"))
    {
        return;
    }

    const dictionary& modifiersDict = profileDict.subDict("modifiers");
    const wordList keys = modifiersDict.toc();

    modifiers_.setSize(keys.size());

    forAll(keys, i)
    {
        const dictionary& modifierDict = modifiersDict.subDict(keys[i]);
        word type = keys[i];

        if (modifierDict.found("type"))
        {
            modifierDict.lookup("type") >> type;
        }

        autoPtr<HeatSourceProfileModifier> modifier
        (
            makeHeatSourceProfileModifier
            (
                type,
                modifierDict,
                laserDict,
                defaultScale,
                cartesianHalfWidths
            )
        );

        modifiers_.set(i, modifier.ptr());
    }
}


scalar HeatSourceProfile::applyModifiers
(
    const scalar baseIntensity,
    const vector& localPosition,
    const scalar areaScale
) const
{
    scalar intensity = max(baseIntensity, scalar(0));

    forAll(modifiers_, i)
    {
        intensity *= max(modifiers_[i].evaluate(localPosition, areaScale), scalar(0));
    }

    return max(intensity, scalar(0));
}


scalar HeatSourceProfile::evaluate
(
    const vector& localPosition,
    const scalar areaScale
) const
{
    return applyModifiers(evaluateBase(localPosition, areaScale), localPosition, areaScale);
}


void HeatSourceProfile::update(const scalar currentTime, const scalar deltaT)
{
    forAll(modifiers_, i)
    {
        modifiers_[i].update(currentTime, deltaT);
    }
}


SuperGaussianProfile::SuperGaussianProfile
(
    const dictionary& dict,
    const dictionary& laserDict
)
:
    HeatSourceProfile(dict),
    order_
    (
        dict.lookupOrDefault<scalar>
        (
            "order",
            dict.lookupOrDefault<scalar>("order", scalar(2))
        )
    ),
    cutoff_(dict.lookupOrDefault<scalar>("cutoff", scalar(1e-6))),
    supportRadius_(SMALL),
    centerShift_
    (
        dict.lookupOrDefault<vector>
        (
            "centerShift",
            dict.lookupOrDefault<vector>("shift", vector::zero)
        )
    ),
    qXX_(0),
    qXZ_(0),
    qZZ_(0)
{
    scalar defaultRadius = dict.lookupOrDefault<scalar>("radius", scalar(-1));

    if (defaultRadius <= SMALL)
    {
        defaultRadius =
            dict.lookupOrDefault<scalar>
            (
                "laserRadius",
                laserDict.lookupOrDefault<scalar>("laserRadius", scalar(-1))
            );
    }

    if (defaultRadius <= SMALL)
    {
        defaultRadius = dict.lookupOrDefault<scalar>("HS_a", scalar(-1));
    }

    scalar radiusX =
        dict.lookupOrDefault<scalar>
        (
            "radiusX",
            dict.lookupOrDefault<scalar>("radius", defaultRadius)
        );
    scalar radiusZ =
        dict.lookupOrDefault<scalar>
        (
            "radiusZ",
            dict.lookupOrDefault<scalar>("radius", defaultRadius)
        );
    const scalar rotation =
        degToRad(dict.lookupOrDefault<scalar>("rotation", scalar(0)));

    if (order_ <= SMALL)
    {
        FatalErrorInFunction
            << "SuperGaussianProfile expects 'order' > 0"
            << nl << dict << exit(FatalError);
    }

    cutoff_ = min(max(cutoff_, scalar(1e-12)), scalar(0.999999));

    if (dict.found("quadraticForm"))
    {
        readSymmetric2x2(dict.subDict("quadraticForm"), qXX_, qXZ_, qZZ_);
    }
    else if (dict.found("inverseCovariance"))
    {
        readSymmetric2x2(dict.subDict("inverseCovariance"), qXX_, qXZ_, qZZ_);
    }
    else if (dict.found("covarianceMatrix") || dict.found("varianceMatrix"))
    {
        const dictionary& covDict =
            dict.found("covarianceMatrix")
          ? dict.subDict("covarianceMatrix")
          : dict.subDict("varianceMatrix");

        scalar covXX = 0.0;
        scalar covXZ = 0.0;
        scalar covZZ = 0.0;
        readSymmetric2x2(covDict, covXX, covXZ, covZZ);
        invertSymmetric2x2(covXX, covXZ, covZZ, qXX_, qXZ_, qZZ_);
    }
    else
    {
        scalar majorRadius =
            dict.lookupOrDefault<scalar>("majorRadius", radiusX);
        scalar minorRadius =
            dict.lookupOrDefault<scalar>("minorRadius", radiusZ);

        if
        (
            dict.found("majorVariance")
         || dict.found("minorVariance")
        )
        {
            majorRadius =
                Foam::sqrt
                (
                    max
                    (
                        dict.lookupOrDefault<scalar>("majorVariance", sqr(majorRadius)),
                        SMALL
                    )
                );
            minorRadius =
                Foam::sqrt
                (
                    max
                    (
                        dict.lookupOrDefault<scalar>("minorVariance", sqr(minorRadius)),
                        SMALL
                    )
                );
        }

        if (majorRadius <= SMALL || minorRadius <= SMALL)
        {
            FatalErrorInFunction
                << "SuperGaussianProfile requires positive radius or variance data."
                << nl << dict << exit(FatalError);
        }

        const scalar c = Foam::cos(rotation);
        const scalar s = Foam::sin(rotation);
        const scalar invMajor2 = scalar(1)/sqr(majorRadius);
        const scalar invMinor2 = scalar(1)/sqr(minorRadius);

        qXX_ = c*c*invMajor2 + s*s*invMinor2;
        qXZ_ = c*s*(invMajor2 - invMinor2);
        qZZ_ = s*s*invMajor2 + c*c*invMinor2;
    }

    const scalar lambdaMin = smallestEigenvalue2x2(qXX_, qXZ_, qZZ_);

    if (lambdaMin <= SMALL)
    {
        FatalErrorInFunction
            << "SuperGaussianProfile expects a positive-definite quadratic form."
            << nl << dict << exit(FatalError);
    }

    const scalar alpha = Foam::pow(-Foam::log(cutoff_), scalar(1)/order_);
    const scalar effectiveRadius = scalar(1)/Foam::sqrt(lambdaMin);

    supportRadius_ =
        alpha*effectiveRadius
      + Foam::sqrt(sqr(centerShift_.x()) + sqr(centerShift_.z()));

    buildModifiers
    (
        dict,
        laserDict,
        defaultProfileScale(dict, laserDict, supportRadius_),
        cartesianSupportHalfWidths()
    );
}


AnnularProfile::AnnularProfile
(
    const dictionary& dict,
    const dictionary& laserDict
)
:
    HeatSourceProfile(dict),
    innerRadius_(scalar(0)),
    outerRadius_(scalar(0)),
    falloff_(scalar(0)),
    rotation_(degToRad(dict.lookupOrDefault<scalar>("rotation", scalar(0)))),
    peakAtInner_(Switch(dict.lookupOrDefault<bool>("peakAtInner", false))),
    supportRadius_(SMALL)
{
    dict.lookup("innerRadius") >> innerRadius_;
    dict.lookup("outerRadius") >> outerRadius_;

    if (innerRadius_ < SMALL || outerRadius_ <= innerRadius_)
    {
        FatalErrorInFunction
            << "AnnularProfile expects 0 <= innerRadius < outerRadius"
            << nl << dict << exit(FatalError);
    }

    falloff_ = dict.lookupOrDefault<scalar>
    (
        "falloff",
        scalar(0.1)*(outerRadius_ - innerRadius_)
    );

    falloff_ = max(falloff_, SMALL);
    supportRadius_ = outerRadius_ + scalar(5)*falloff_;

    buildModifiers
    (
        dict,
        laserDict,
        supportRadius_,
        cartesianSupportHalfWidths()
    );
}


CsvProfile::CsvProfile
(
    const dictionary& dict,
    const dictionary& laserDict
)
:
    HeatSourceProfile(dict),
    file_(dict.get<fileName>("file").expand()),
    nx_(0),
    nz_(0),
    dx_(0),
    dz_(0),
    halfX_(0),
    halfZ_(0),
    values_(0)
{
    dict.lookup("nx") >> nx_;
    nz_ = dict.lookupOrDefault<label>("nz", nx_);
    dict.lookup("dx") >> dx_;
    dz_ = dict.lookupOrDefault<scalar>("dz", dict.lookupOrDefault<scalar>("dy", dx_));

    if (nx_ <= 1 || nz_ <= 1)
    {
        FatalErrorInFunction
            << "CsvProfile requires nx > 1 and nz > 1"
            << nl << dict << exit(FatalError);
    }

    halfX_ = scalar(0.5)*dict.lookupOrDefault<scalar>
    (
        "xLength",
        dx_*scalar(nx_ - 1)
    );
    halfZ_ = scalar(0.5)*dict.lookupOrDefault<scalar>
    (
        "zLength",
        dz_*scalar(nz_ - 1)
    );

    values_.setSize(nx_*nz_, scalar(0));

    IFstream is(file_);
    if (!is.good())
    {
        FatalErrorInFunction
            << "Unable to open CSV intensity file '" << file_ << "'"
            << nl << dict << exit(FatalError);
    }

    forAll(values_, i)
    {
        scalar val = 0.0;
        is >> val;

        if (!is.good())
        {
            FatalErrorInFunction
                << "Failed to read entry " << i << " from intensity file '"
                << file_ << "'"
                << nl << dict << exit(FatalError);
        }

        values_[i] = max(val, scalar(0));
    }

    buildModifiers
    (
        dict,
        laserDict,
        supportRadius(),
        cartesianSupportHalfWidths()
    );
}


BesselProfile::BesselProfile
(
    const dictionary& dict,
    const dictionary& laserDict
)
:
    HeatSourceProfile(dict),
    radius_(scalar(0)),
    firstZeroRadius_(scalar(-1)),
    order_(dict.lookupOrDefault<label>("order", 0)),
    cutoff_(dict.lookupOrDefault<scalar>("cutoff", scalar(1e-3))),
    supportRadius_(SMALL)
{
    if (order_ < 0)
    {
        FatalErrorInFunction
            << "BesselProfile expects 'order' >= 0"
            << nl << dict << exit(FatalError);
    }

    cutoff_ = min(max(cutoff_, scalar(1e-12)), scalar(0.999999));

    scalar estimatedFirstZero = estimateBesselZero(order_, 1);

    firstZeroRadius_ =
        dict.lookupOrDefault<scalar>("firstZeroRadius", scalar(-1));

    if (firstZeroRadius_ > SMALL)
    {
        radius_ = firstZeroRadius_/estimatedFirstZero;
    }
    else
    {
        scalar defaultRadius = dict.lookupOrDefault<scalar>("radius", scalar(-1));

        if (defaultRadius <= SMALL)
        {
            defaultRadius =
                dict.lookupOrDefault<scalar>
                (
                    "laserRadius",
                    laserDict.lookupOrDefault<scalar>("laserRadius", scalar(-1))
                );
        }

        radius_ = defaultRadius;
    }

    if (radius_ <= SMALL)
    {
        FatalErrorInFunction
            << "BesselProfile requires positive radius or firstZeroRadius value"
            << nl << dict << exit(FatalError);
    }

    const scalar estimatedSupportRadius =
        radius_*estimatedFirstZero*(scalar(1) + Foam::log(scalar(1)/cutoff_));

    supportRadius_ =
        dict.lookupOrDefault<scalar>("supportRadius", estimatedSupportRadius);

    if (supportRadius_ <= SMALL)
    {
        FatalErrorInFunction
            << "BesselProfile expects positive supportRadius when specified"
            << nl << dict << exit(FatalError);
    }

    supportRadius_ = max(supportRadius_, estimatedFirstZero*radius_);

    buildModifiers
    (
        dict,
        laserDict,
        supportRadius_,
        cartesianSupportHalfWidths()
    );
}


#if LASERHEATSOURCE_HAS_BEAMSHAPER
BeamShapeProfile::BeamShapeProfile
(
    const dictionary& dict,
    const dictionary& laserDict,
    const fvMesh& mesh,
    const dictionary& laserProperties
)
:
    HeatSourceProfile(dict),
    beamShape_
    (
        new BeamShape
        (
            dict.lookupOrDefault<fileName>
            (
                "shapeFile",
                fileName("beamShape.inp")
            ),
            dict.lookupOrDefault<fileName>
            (
                "scheduleFile",
                fileName("beamShapeDict")
            ),
            mesh,
            laserProperties
        )
    ),
    halfLength_(scalar(0.5)*beamShape_().imageLength_.value()),
    cacheShapePerStep_(dict.lookupOrDefault<bool>("cachePerStep", true))
{
    buildModifiers
    (
        dict,
        laserDict,
        halfLength_,
        cartesianSupportHalfWidths()
    );
}


scalar SuperGaussianProfile::evaluateBase
(
    const vector& localPosition,
    const scalar
) const
{
    const vector shifted = localPosition - centerShift_;
    const scalar xLocal = shifted.x();
    const scalar zLocal = shifted.z();

    const scalar quadratic =
        qXX_*sqr(xLocal) + scalar(2)*qXZ_*xLocal*zLocal + qZZ_*sqr(zLocal);

    const scalar radius = Foam::sqrt(max(quadratic, scalar(0)));

    return max(Foam::exp(-Foam::pow(radius, order_)), scalar(0));
}


scalar AnnularProfile::evaluateBase
(
    const vector& localPosition,
    const scalar
) const
{
    const scalar c = Foam::cos(rotation_);
    const scalar s = Foam::sin(rotation_);

    const scalar xLocal = c*localPosition.x() + s*localPosition.z();
    const scalar zLocal = -s*localPosition.x() + c*localPosition.z();

    const scalar r = Foam::sqrt(sqr(xLocal) + sqr(zLocal));

    const scalar innerTransition =
        scalar(0.5)*(scalar(1) + Foam::tanh((r - innerRadius_)/falloff_));
    const scalar outerTransition =
        scalar(0.5)*(scalar(1) - Foam::tanh((r - outerRadius_)/falloff_));

    scalar intensity = innerTransition*outerTransition;

    if (peakAtInner_)
    {
        intensity *= innerRadius_/max(r, SMALL);
    }

    return max(intensity, scalar(0));
}


scalar CsvProfile::evaluateBase
(
    const vector& localPosition,
    const scalar
) const
{
    const scalar x = localPosition.x();
    const scalar z = localPosition.z();

    if (mag(x) > halfX_ + dx_ || mag(z) > halfZ_ + dz_)
    {
        return 0.0;
    }

    scalar fx = (x + halfX_)/dx_;
    scalar fz = (z + halfZ_)/dz_;

    fx = min(max(fx, scalar(0)), scalar(nx_ - 1));
    fz = min(max(fz, scalar(0)), scalar(nz_ - 1));

    const label ix0 = min(label(fx), nx_ - 1);
    const label iz0 = min(label(fz), nz_ - 1);
    const label ix1 = min(ix0 + 1, nx_ - 1);
    const label iz1 = min(iz0 + 1, nz_ - 1);
    const scalar tx = fx - scalar(ix0);
    const scalar tz = fz - scalar(iz0);

    const scalar v00 = values_[iz0*nx_ + ix0];
    const scalar v10 = values_[iz0*nx_ + ix1];
    const scalar v01 = values_[iz1*nx_ + ix0];
    const scalar v11 = values_[iz1*nx_ + ix1];
    const scalar interpX0 = v00 + tx*(v10 - v00);
    const scalar interpX1 = v01 + tx*(v11 - v01);

    return max(interpX0 + tz*(interpX1 - interpX0), scalar(0));
}


scalar BesselProfile::evaluateBase
(
    const vector& localPosition,
    const scalar
) const
{
    const scalar r =
        Foam::sqrt(sqr(localPosition.x()) + sqr(localPosition.z()));

    if (r > supportRadius_)
    {
        return 0.0;
    }

    const scalar x = r/radius_;
    const scalar intensity = besselFunction(order_, x);

    return max(sqr(intensity), scalar(0));
}


scalar BeamShapeProfile::evaluateBase
(
    const vector& localPosition,
    const scalar localScale
) const
{
    if (!beamShape_.valid())
    {
        return 0.0;
    }

    return max
    (
        beamShape_->interpolateValue
        (
            localPosition.x(),
            localPosition.z(),
            max(localScale, SMALL),
            550
        ),
        scalar(0)
    );
}


void BeamShapeProfile::update(const scalar currentTime, const scalar deltaT)
{
    HeatSourceProfile::update(currentTime, deltaT);

    if (beamShape_.valid())
    {
        beamShape_->getBeamShape(currentTime, deltaT);
    }
}
#endif


HeatSourceProfileSet::HeatSourceProfileSet()
:
    profiles_(0)
{}


HeatSourceProfileSet::HeatSourceProfileSet(const label nProfiles)
:
    profiles_(nProfiles)
{}


void HeatSourceProfileSet::setProfile(const label i, HeatSourceProfile* profile)
{
    profiles_.set(i, profile);
}


void HeatSourceProfileSet::update(const scalar currentTime, const scalar deltaT)
{
    forAll(profiles_, i)
    {
        profiles_[i].update(currentTime, deltaT);
    }
}


autoPtr<HeatSourceProfile> makeHeatSourceProfile
(
    const word& rawType,
    const dictionary& dict,
    const dictionary& laserDict,
    const fvMesh& mesh,
    const dictionary& laserProperties
)
{
    word type(rawType);
    type = word(stringOps::lower(type));

    if (type == "gaussian" || type == "supergaussian")
    {
        return autoPtr<HeatSourceProfile>(new SuperGaussianProfile(dict, laserDict));
    }
    else if (type == "annular" || type == "annulus")
    {
        return autoPtr<HeatSourceProfile>(new AnnularProfile(dict, laserDict));
    }
    else if (type == "csv")
    {
        return autoPtr<HeatSourceProfile>(new CsvProfile(dict, laserDict));
    }
    else if (type == "bessel")
    {
        return autoPtr<HeatSourceProfile>(new BesselProfile(dict, laserDict));
    }
    else if (type == "beamshape" || type == "beamshaper")
    {
#if LASERHEATSOURCE_HAS_BEAMSHAPER
        return autoPtr<HeatSourceProfile>
        (
            new BeamShapeProfile(dict, laserDict, mesh, laserProperties)
        );
#else
        FatalErrorInFunction
            << "The 'beamShape' profile type is not available in this public "
            << "build of laserHeatSource."
            << nl
            << "Other profile types remain available, but beam-shaper support "
            << "has been excluded from this source tree."
            << nl << dict << exit(FatalError);

        return autoPtr<HeatSourceProfile>();
#endif
    }

    FatalErrorInFunction
        << "Unknown heat source profile type '" << rawType << "'"
        << nl << dict << exit(FatalError);

    return autoPtr<HeatSourceProfile>();
}


autoPtr<HeatSourceProfileSet> buildHeatSourceProfileSet
(
    const dictionary& laserDict,
    const fvMesh& mesh,
    const dictionary& laserProperties
)
{
    if (laserDict.found("profile"))
    {
        FatalErrorInFunction
            << "Use a 'profiles' dictionary instead of 'profile'."
            << nl << laserDict << exit(FatalError);
    }

    if (!laserDict.found("profiles"))
    {
        return autoPtr<HeatSourceProfileSet>(new HeatSourceProfileSet());
    }

    const dictionary& profilesDict = laserDict.subDict("profiles");
    const wordList keys = profilesDict.toc();

    if (keys.empty())
    {
        FatalErrorInFunction
            << "The 'profiles' dictionary must contain at least one profile."
            << nl << laserDict << exit(FatalError);
    }

    autoPtr<HeatSourceProfileSet> setPtr
    (
        new HeatSourceProfileSet(keys.size())
    );

    forAll(keys, i)
    {
        const dictionary& profileDict = profilesDict.subDict(keys[i]);
        word type = keys[i];

        if (profileDict.found("type"))
        {
            profileDict.lookup("type") >> type;
        }

        autoPtr<HeatSourceProfile> profile
        (
            makeHeatSourceProfile(type, profileDict, laserDict, mesh, laserProperties)
        );

        setPtr->setProfile(i, profile.ptr());
    }

    return setPtr;
}

} // End namespace Foam

// ************************************************************************* //
