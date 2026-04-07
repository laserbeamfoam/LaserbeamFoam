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

Class
    Foam::HeatSourceProfile and derived classes

Description
    Hierarchy of beam intensity profiles that can be combined to describe
    arbitrary laser heat-source distributions.

\*---------------------------------------------------------------------------*/

#ifndef HeatSourceProfiles_H
#define HeatSourceProfiles_H

#include "List.H"
#include "PtrList.H"
#include "Random.H"
#include "Switch.H"
#include "autoPtr.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "scalarField.H"
#include "vector.H"
#if defined(__has_include)
#if __has_include("beamShaper.H")
#include "beamShaper.H"
#define LASERHEATSOURCE_HAS_BEAMSHAPER 1
#else
#define LASERHEATSOURCE_HAS_BEAMSHAPER 0
#endif
#else
#include "beamShaper.H"
#define LASERHEATSOURCE_HAS_BEAMSHAPER 1
#endif

namespace Foam
{

enum class SamplingStrategy
{
    RADIAL_POLAR,
    CARTESIAN
};


class HeatSourceProfileSampling
{
    SamplingStrategy type_;
    scalar radius_;
    label nRadial_;
    label nAngular_;
    scalar halfWidthX_;
    scalar halfWidthZ_;
    label nCartesianX_;
    label nCartesianZ_;
    scalar minRayPowerFraction_;

public:

    HeatSourceProfileSampling();

    explicit HeatSourceProfileSampling(const dictionary& profileDict);

    SamplingStrategy type() const
    {
        return type_;
    }

    scalar radius() const
    {
        return radius_;
    }

    label nRadial() const
    {
        return nRadial_;
    }

    label nAngular() const
    {
        return nAngular_;
    }

    scalar halfWidthX() const
    {
        return halfWidthX_;
    }

    scalar halfWidthZ() const
    {
        return halfWidthZ_;
    }

    label nCartesianX() const
    {
        return nCartesianX_;
    }

    label nCartesianZ() const
    {
        return nCartesianZ_;
    }

    scalar minRayPowerFraction() const
    {
        return minRayPowerFraction_;
    }
};


class HeatSourceProfileModifier
{
public:

    virtual ~HeatSourceProfileModifier() = default;

    virtual scalar evaluate
    (
        const vector& localPosition,
        const scalar areaScale
    ) const = 0;

    virtual void update(const scalar currentTime, const scalar deltaT)
    {}
};


class PlaneBiasModifier
:
    public HeatSourceProfileModifier
{
    vector localDirection_;
    vector localOrigin_;
    scalar amplitude_;
    scalar offset_;
    scalar referenceScale_;
    scalar minFactor_;

public:

    PlaneBiasModifier
    (
        const dictionary& dict,
        const dictionary& laserDict,
        const scalar defaultScale
    );

    virtual scalar evaluate
    (
        const vector& localPosition,
        const scalar areaScale
    ) const;
};


struct DrumNoiseMode
{
    label azimuthalOrder;
    label radialOrder;
    scalar amplitude;
    scalar phase;
    scalar frequency;
    scalar currentAmplitude;
    scalar currentPhase;
    scalar normalisation;
    scalar meanValue;

    DrumNoiseMode();
};


struct CartesianNoiseMode
{
    label xOrder;
    label zOrder;
    scalar amplitude;
    scalar phase;
    scalar frequency;
    scalar currentAmplitude;
    scalar currentPhase;

    CartesianNoiseMode();
};


class DrumNoiseModifier
:
    public HeatSourceProfileModifier
{
    List<DrumNoiseMode> modes_;
    scalar radius_;
    scalar minFactor_;
    Switch randomPhase_;
    Switch randomAmplitude_;
    scalar amplitudeJitter_;
    word evolution_;
    mutable scalar lastUpdateTime_;
    Random rndGen_;

    scalar modeSignal
    (
        const DrumNoiseMode& mode,
        const scalar r,
        const scalar theta
    ) const;

public:

    DrumNoiseModifier
    (
        const dictionary& dict,
        const scalar defaultRadius
    );

    virtual scalar evaluate
    (
        const vector& localPosition,
        const scalar areaScale
    ) const;

    virtual void update(const scalar currentTime, const scalar deltaT);
};


class CartesianModeNoiseModifier
:
    public HeatSourceProfileModifier
{
    List<CartesianNoiseMode> modes_;
    scalar halfWidth_;
    scalar halfHeight_;
    scalar minFactor_;
    Switch randomPhase_;
    Switch randomAmplitude_;
    scalar amplitudeJitter_;
    word evolution_;
    mutable scalar lastUpdateTime_;
    Random rndGen_;

    scalar modeSignal
    (
        const CartesianNoiseMode& mode,
        const scalar x,
        const scalar z
    ) const;

public:

    CartesianModeNoiseModifier
    (
        const dictionary& dict,
        const vector& defaultHalfWidths
    );

    virtual scalar evaluate
    (
        const vector& localPosition,
        const scalar areaScale
    ) const;

    virtual void update(const scalar currentTime, const scalar deltaT);
};


class GridNoiseModifier
:
    public HeatSourceProfileModifier
{
    label nx_;
    label nz_;
    scalar halfWidth_;
    scalar halfHeight_;
    scalar amplitude_;
    scalar minFactor_;
    scalar correlationLength_;
    label smoothPasses_;
    scalar frequency_;
    word evolution_;
    mutable scalar currentBlend_;
    mutable scalar anchorTime_;
    mutable scalar nextAnchorTime_;
    mutable scalar lastUpdateTime_;
    mutable scalarField currentMask_;
    mutable scalarField nextMask_;
    Random rndGen_;

    void generateMask(scalarField& mask);

    scalar sampleMask
    (
        const scalarField& mask,
        const scalar x,
        const scalar z
    ) const;

public:

    GridNoiseModifier
    (
        const dictionary& dict,
        const vector& defaultHalfWidths
    );

    virtual scalar evaluate
    (
        const vector& localPosition,
        const scalar areaScale
    ) const;

    virtual void update(const scalar currentTime, const scalar deltaT);
};


class HeatSourceProfile
{
protected:

    scalar weight_;
    HeatSourceProfileSampling sampling_;
    PtrList<HeatSourceProfileModifier> modifiers_;

    void buildModifiers
    (
        const dictionary& profileDict,
        const dictionary& laserDict,
        const scalar defaultScale,
        const vector& cartesianHalfWidths
    );

    scalar applyModifiers
    (
        const scalar baseIntensity,
        const vector& localPosition,
        const scalar areaScale
    ) const;

public:

    explicit HeatSourceProfile(const dictionary& dict);

    virtual ~HeatSourceProfile() = default;

    virtual scalar evaluateBase
    (
        const vector& localPosition,
        const scalar areaScale
    ) const = 0;

    virtual scalar evaluate
    (
        const vector& localPosition,
        const scalar areaScale
    ) const;

    virtual scalar supportRadius() const = 0;

    virtual vector cartesianSupportHalfWidths() const
    {
        const scalar r = supportRadius();
        return vector(r, 0, r);
    }

    virtual void update(const scalar currentTime, const scalar deltaT);

    scalar weight() const
    {
        return weight_;
    }

    const HeatSourceProfileSampling& sampling() const
    {
        return sampling_;
    }
};


class SuperGaussianProfile
:
    public HeatSourceProfile
{
    scalar order_;
    scalar cutoff_;
    scalar supportRadius_;
    vector centerShift_;
    scalar qXX_;
    scalar qXZ_;
    scalar qZZ_;

public:

    SuperGaussianProfile
    (
        const dictionary& dict,
        const dictionary& laserDict
    );

    virtual scalar evaluateBase
    (
        const vector& localPosition,
        const scalar areaScale
    ) const;

    virtual scalar supportRadius() const
    {
        return supportRadius_;
    }

};


class AnnularProfile
:
    public HeatSourceProfile
{
    scalar innerRadius_;
    scalar outerRadius_;
    scalar falloff_;
    scalar rotation_;
    Switch peakAtInner_;
    scalar supportRadius_;

public:

    AnnularProfile
    (
        const dictionary& dict,
        const dictionary& laserDict
    );

    virtual scalar evaluateBase
    (
        const vector& localPosition,
        const scalar areaScale
    ) const;

    virtual scalar supportRadius() const
    {
        return supportRadius_;
    }

};


class CsvProfile
:
    public HeatSourceProfile
{
    fileName file_;
    label nx_;
    label nz_;
    scalar dx_;
    scalar dz_;
    scalar halfX_;
    scalar halfZ_;
    scalarField values_;

public:

    CsvProfile
    (
        const dictionary& dict,
        const dictionary& laserDict
    );

    virtual scalar evaluateBase
    (
        const vector& localPosition,
        const scalar areaScale
    ) const;

    virtual scalar supportRadius() const
    {
        return Foam::sqrt(sqr(halfX_) + sqr(halfZ_));
    }

    virtual vector cartesianSupportHalfWidths() const
    {
        return vector(halfX_, 0, halfZ_);
    }
};


class BesselProfile
:
    public HeatSourceProfile
{
    scalar radius_;
    scalar firstZeroRadius_;
    label order_;
    scalar cutoff_;
    scalar supportRadius_;

public:

    BesselProfile
    (
        const dictionary& dict,
        const dictionary& laserDict
    );

    virtual scalar evaluateBase
    (
        const vector& localPosition,
        const scalar areaScale
    ) const;

    virtual scalar supportRadius() const
    {
        return supportRadius_;
    }

};


#if LASERHEATSOURCE_HAS_BEAMSHAPER
class BeamShapeProfile
:
    public HeatSourceProfile
{
    mutable autoPtr<BeamShape> beamShape_;
    scalar halfLength_;
    Switch cacheShapePerStep_;

public:

    BeamShapeProfile
    (
        const dictionary& dict,
        const dictionary& laserDict,
        const fvMesh& mesh,
        const dictionary& laserProperties
    );

    virtual scalar evaluateBase
    (
        const vector& localPosition,
        const scalar areaScale
    ) const;

    virtual scalar supportRadius() const
    {
        return halfLength_;
    }

    virtual vector cartesianSupportHalfWidths() const
    {
        return vector(halfLength_, 0, halfLength_);
    }

    virtual void update(const scalar currentTime, const scalar deltaT);
};
#endif


class HeatSourceProfileSet
{
    PtrList<HeatSourceProfile> profiles_;

public:

    HeatSourceProfileSet();

    explicit HeatSourceProfileSet(const label nProfiles);

    void setProfile(const label i, HeatSourceProfile* profile);

    bool empty() const
    {
        return profiles_.empty();
    }

    label size() const
    {
        return profiles_.size();
    }

    void update(const scalar currentTime, const scalar deltaT);

    const PtrList<HeatSourceProfile>& profiles() const
    {
        return profiles_;
    }
};


autoPtr<HeatSourceProfile> makeHeatSourceProfile
(
    const word& rawType,
    const dictionary& dict,
    const dictionary& laserDict,
    const fvMesh& mesh,
    const dictionary& laserProperties
);


autoPtr<HeatSourceProfileSet> buildHeatSourceProfileSet
(
    const dictionary& laserDict,
    const fvMesh& mesh,
    const dictionary& laserProperties
);

} // End namespace Foam

#endif

// ************************************************************************* //
