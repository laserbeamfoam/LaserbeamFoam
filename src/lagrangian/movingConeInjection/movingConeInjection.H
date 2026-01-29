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

#ifndef movingConeInjection_H
#define movingConeInjection_H

#include "InjectionModel.H"
#include "distributionModel.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class CloudType>
class movingConeInjection
:
    public InjectionModel<CloudType>
{
    // Private data

        //- Initial/reference position
        vector position_;

        //- Initial/reference direction
        vector direction_;

        //- Trajectory functions for position (optional)
        autoPtr<Function1<scalar>> positionX_;
        autoPtr<Function1<scalar>> positionY_;
        autoPtr<Function1<scalar>> positionZ_;

        //- Cell label
        label injectorCell_;

        //- TetFace label
        label injectorTetFace_;

        //- TetPt label
        label injectorTetPt_;

        //- Injection duration [s]
        scalar duration_;

        //- Number of parcels per second (optional)
        const scalar parcelsPerSecond_;

        //- Number of parcels per injector (optional)
        const scalar parcelsPerInjector_;

        //- Flow rate profile
        autoPtr<Function1<scalar>> flowRateProfile_;

        //- Parcel velocity magnitude [m/s]
        autoPtr<Function1<scalar>> Umag_;

        //- Inner half-cone angle [deg]
        autoPtr<Function1<scalar>> thetaInner_;

        //- Outer half-cone angle [deg]
        autoPtr<Function1<scalar>> thetaOuter_;

        //- Parcel size distribution model
        autoPtr<distributionModel> sizeDistribution_;

        //- Number of parcels injected
        mutable label nInjected_;

        //- First tangential vector
        vector tanVec1_;

        //- Second tangential vector
        vector tanVec2_;


    // Private Member Functions

        //- Get current injection position
        vector currentPosition(const scalar t) const;


public:

    //- Runtime type information
    TypeName("movingConeInjection");


    // Constructors

        //- Construct from dictionary
        movingConeInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        movingConeInjection(const movingConeInjection<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new movingConeInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~movingConeInjection() = default;


    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Set the injection position and owner cell
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFacei,
                label& tetPti
            );

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "movingConeInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //