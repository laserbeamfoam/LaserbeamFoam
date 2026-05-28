/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
\*---------------------------------------------------------------------------*/

#include "movingConeInjection.H"
#include "Function1.H"
#include "mathematicalConstants.H"
#include "unitConversion.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class CloudType>
Foam::vector Foam::movingConeInjection<CloudType>::currentPosition
(
    const scalar t
) const
{
    vector pos = position_;
    
    if (positionX_.valid())
    {
        pos.x() = positionX_->value(t);
    }
    if (positionY_.valid())
    {
        pos.y() = positionY_->value(t);
    }
    if (positionZ_.valid())
    {
        pos.z() = positionZ_->value(t);
    }
    
    return pos;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::movingConeInjection<CloudType>::movingConeInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName, typeName),
    position_(this->coeffDict().lookup("position")),
    direction_(this->coeffDict().lookup("direction")),
    positionX_(nullptr),
    positionY_(nullptr),
    positionZ_(nullptr),
    injectorCell_(-1),
    injectorTetFace_(-1),
    injectorTetPt_(-1),
    duration_(this->coeffDict().getScalar("duration")),
    parcelsPerSecond_
    (
        this->coeffDict().lookupOrDefault("parcelsPerSecond", scalar(-1.0))
    ),
    parcelsPerInjector_
    (
        this->coeffDict().lookupOrDefault("parcelsPerInjector", scalar(-1.0))
    ),
    flowRateProfile_
    (
        Function1<scalar>::New
        (
            "flowRateProfile",
            this->coeffDict(),
            &owner.mesh()
        )
    ),
    Umag_
    (
        Function1<scalar>::New
        (
            "Umag",
            this->coeffDict(),
            &owner.mesh()
        )
    ),
    thetaInner_
    (
        Function1<scalar>::New
        (
            "thetaInner",
            this->coeffDict(),
            &owner.mesh()
        )
    ),
    thetaOuter_
    (
        Function1<scalar>::New
        (
            "thetaOuter",
            this->coeffDict(),
            &owner.mesh()
        )
    ),
    sizeDistribution_
    (
        distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"), owner.rndGen()
        )
    ),
    nInjected_(Pstream::master() ? this->parcelsAddedTotal() : 0),
    tanVec1_(),
    tanVec2_()
{
    updateMesh();

    // Convert from user time to reduce the number of time conversion calls
    const Time& time = owner.db().time();
    duration_ = time.userTimeToTime(duration_);
    flowRateProfile_->userTimeToTime(time);
    Umag_->userTimeToTime(time);
    thetaInner_->userTimeToTime(time);
    thetaOuter_->userTimeToTime(time);

    // Read trajectory functions if present
    if (this->coeffDict().found("trajectory"))
    {
        const dictionary& trajDict = this->coeffDict().subDict("trajectory");
        
        if (trajDict.found("x"))
        {
            positionX_ = Function1<scalar>::New("x", trajDict, &owner.mesh());
            positionX_->userTimeToTime(time);
        }
        if (trajDict.found("y"))
        {
            positionY_ = Function1<scalar>::New("y", trajDict, &owner.mesh());
            positionY_->userTimeToTime(time);
        }
        if (trajDict.found("z"))
        {
            positionZ_ = Function1<scalar>::New("z", trajDict, &owner.mesh());
            positionZ_->userTimeToTime(time);
        }
    }

    // Determine injection mode
    if (parcelsPerInjector_ > 0)
    {
        Info<< "    Using parcelsPerInjector mode: " 
            << parcelsPerInjector_ << " parcels total" << endl;
    }
    else if (parcelsPerSecond_ > 0)
    {
        Info<< "    Using parcelsPerSecond mode: " 
            << parcelsPerSecond_ << " parcels/second" << endl;
    }
    else
    {
        FatalErrorInFunction
            << "Must specify either parcelsPerInjector or parcelsPerSecond"
            << abort(FatalError);
    }

    // Normalise direction vector and determine tangential vectors
    direction_.normalise();

    vector tangent = Zero;
    scalar magTangent = 0.0;

    Random& rnd = this->owner().rndGen();
    while (magTangent < SMALL)
    {
        vector v = rnd.sample01<vector>();
        tangent = v - (v & direction_)*direction_;
        magTangent = mag(tangent);
    }

    tanVec1_ = tangent/magTangent;
    tanVec2_ = direction_^tanVec1_;

    // Set total volume to inject
    this->volumeTotal_ = flowRateProfile_->integrate(0.0, duration_);
}


template<class CloudType>
Foam::movingConeInjection<CloudType>::movingConeInjection
(
    const movingConeInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    position_(im.position_),
    direction_(im.direction_),
    positionX_(im.positionX_.clone()),
    positionY_(im.positionY_.clone()),
    positionZ_(im.positionZ_.clone()),
    injectorCell_(im.injectorCell_),
    injectorTetFace_(im.injectorTetFace_),
    injectorTetPt_(im.injectorTetPt_),
    duration_(im.duration_),
    parcelsPerSecond_(im.parcelsPerSecond_),
    parcelsPerInjector_(im.parcelsPerInjector_),
    flowRateProfile_(im.flowRateProfile_.clone()),
    Umag_(im.Umag_.clone()),
    thetaInner_(im.thetaInner_.clone()),
    thetaOuter_(im.thetaOuter_.clone()),
    sizeDistribution_(im.sizeDistribution_.clone()),
    nInjected_(im.nInjected_),
    tanVec1_(im.tanVec1_),
    tanVec2_(im.tanVec2_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::movingConeInjection<CloudType>::updateMesh()
{
    // Get current position (at SOI for initial update)
    vector pos = currentPosition(this->SOI_);
    
    // Set/cache the injector cell
    if
    (
        !this->findCellAtPosition
        (
            injectorCell_,
            injectorTetFace_,
            injectorTetPt_,
            pos,
            !this->ignoreOutOfBounds_
        )
    )
    {
        WarningInFunction
            << "Initial injection position " << pos 
            << " is out of bounds" << endl;
    }
}


template<class CloudType>
Foam::scalar Foam::movingConeInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
Foam::label Foam::movingConeInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        // Calculate target volume based on flow rate profile
        const scalar targetVolume = flowRateProfile_->integrate(0, time1);
        const scalar volumeFraction = targetVolume/this->volumeTotal_;
        
        label targetParcels = 0;
        
        if (parcelsPerInjector_ > 0)
        {
            // Mode 1: Total parcels distributed over time (like standard coneInjection)
            targetParcels = ceil(parcelsPerInjector_*volumeFraction);
        }
        else
        {
            // Mode 2: Parcels per second
            targetParcels = ceil(parcelsPerSecond_*duration_*volumeFraction);
        }

        return targetParcels - returnReduce(nInjected_, sumOp<label>());
    }

    return 0;
}


template<class CloudType>
Foam::scalar Foam::movingConeInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return flowRateProfile_->integrate(time0, time1);
    }

    return 0.0;
}


template<class CloudType>
void Foam::movingConeInjection<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar time,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    // Get current position from trajectory
    position = currentPosition(time);
    
    // Find cell at current position
    this->findCellAtPosition
    (
        cellOwner,
        tetFacei,
        tetPti,
        position,
        false
    );
}


template<class CloudType>
void Foam::movingConeInjection<CloudType>::setProperties
(
    const label,
    const label,
    const scalar time,
    typename CloudType::parcelType& parcel
)
{
    Random& rnd = this->owner().rndGen();

    // Set direction vectors
    scalar t = time - this->SOI_;
    scalar ti = thetaInner_->value(t);
    scalar to = thetaOuter_->value(t);
    scalar coneAngle = degToRad(rnd.position<scalar>(ti, to));

    scalar alpha = sin(coneAngle);
    scalar dcorr = cos(coneAngle);
    scalar beta = twoPi*rnd.sample01<scalar>();

    vector normal = alpha*(tanVec1_*cos(beta) + tanVec2_*sin(beta));
    vector dirVec = dcorr*direction_;
    dirVec += normal;
    dirVec.normalise();

    // Set particle velocity
    parcel.U() = Umag_->value(t)*dirVec;

    // Set particle diameter
    parcel.d() = sizeDistribution_().sample();

    // Increment number of particles injected
    nInjected_++;
}


template<class CloudType>
bool Foam::movingConeInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::movingConeInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //