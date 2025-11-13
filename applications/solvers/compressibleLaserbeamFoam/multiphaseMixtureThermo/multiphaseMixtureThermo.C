/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "multiphaseMixtureThermo.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "Time.H"
#include "subCycle.H"
#include "MULES.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcMeshPhi.H"
#include "surfaceInterpolate.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphaseMixtureThermo, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiphaseMixtureThermo::calcAlphas()
{
    scalar level = 0.0;
    alphas_ == 0.0;

    for (const phaseModel& phase : phases_)
    {
        alphas_ += level * phase;
        level += 1.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseMixtureThermo::multiphaseMixtureThermo
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    psiThermo(U.mesh(), word::null),
    phases_(lookup("phases"), phaseModel::iNew(p_, T_)),

    mesh_(U.mesh()),
    U_(U),
    phi_(phi),

    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimTime, Zero)
    ),

    alphas_
    (
        IOobject
        (
            "alphas",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),

    sigmas_(lookup("sigmas")),
    dimSigma_(1, 0, -2, 0, 0),
    dsigmadT_(lookup("dsigmadT")),
    dimdSigmadT_(1, 0, -2, -1, 0),
    boils_(lookup("boils")),
    dimBoil_(0, 0, 0, 1, 0),
    LatentHeatGass_(lookup("LatentHeatGas")),
    dimLatentHeatGas_(0, 2, -2, 0, 0),
    cAlphas_(lookup("interfaceCompression")),
    dAlphas_(lookup("interfaceDiffusion")),
    dimdiff_(0, 2, -1, 0, 0),
    deltaN_
    (
        "deltaN",
        1e-8/cbrt(average(mesh_.V()))
    )
{
    rhoPhi_.setOriented();
    calcAlphas();
    alphas_.write();
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::multiphaseMixtureThermo::correct()
{
    for (phaseModel& phase : phases_)
    {
        phase.correct();
    }

    auto phasei = phases_.cbegin();

    psi_ = phasei()*phasei().thermo().psi();
    mu_ = phasei()*phasei().thermo().mu();
    alpha_ = phasei()*phasei().thermo().alpha();

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        psi_ += phasei()*phasei().thermo().psi();
        mu_ += phasei()*phasei().thermo().mu();
        alpha_ += phasei()*phasei().thermo().alpha();
    }
}


void Foam::multiphaseMixtureThermo::correctRho(const volScalarField& dp)
{
    for (phaseModel& phase : phases_)
    {
        phase.thermo().rho() += phase.thermo().psi()*dp;
    }
}


Foam::word Foam::multiphaseMixtureThermo::thermoName() const
{
    auto phasei = phases_.cbegin();

    word name = phasei().thermo().thermoName();

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        name += ',' + phasei().thermo().thermoName();
    }

    return name;
}


bool Foam::multiphaseMixtureThermo::incompressible() const
{
    for (const phaseModel& phase : phases_)
    {
        if (!phase.thermo().incompressible())
        {
            return false;
        }
    }

    return true;
}


bool Foam::multiphaseMixtureThermo::isochoric() const
{
    for (const phaseModel& phase : phases_)
    {
        if (!phase.thermo().isochoric())
        {
            return false;
        }
    }

    return true;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> the(phasei()*phasei().thermo().he(p, T));

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        the.ref() += phasei()*phasei().thermo().he(p, T);
    }

    return the;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> the
    (
        scalarField(phasei(), cells)*phasei().thermo().he(p, T, cells)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        the.ref() +=
            scalarField(phasei(), cells)*phasei().thermo().he(p, T, cells);
    }

    return the;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::he
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> the
    (
        phasei().boundaryField()[patchi]*phasei().thermo().he(p, T, patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        the.ref() +=
            phasei().boundaryField()[patchi]*phasei().thermo().he(p, T, patchi);
    }

    return the;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::hc() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> thc(phasei()*phasei().thermo().hc());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        thc.ref() += phasei()*phasei().thermo().hc();
    }

    return thc;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const labelList& cells
) const
{
    NotImplemented;
    return T0;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::THE
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& T0,
    const label patchi
) const
{
    NotImplemented;
    return T0;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::rho() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> trho(phasei()*phasei().thermo().rho());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        trho.ref() += phasei()*phasei().thermo().rho();
    }

    return trho;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::rho
(
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> trho
    (
        phasei().boundaryField()[patchi]*phasei().thermo().rho(patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        trho.ref() +=
            phasei().boundaryField()[patchi]*phasei().thermo().rho(patchi);
    }

    return trho;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::Cp() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> tCp(phasei()*phasei().thermo().Cp());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tCp.ref() += phasei()*phasei().thermo().Cp();
    }

    return tCp;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::Cp
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> tCp
    (
        phasei().boundaryField()[patchi]*phasei().thermo().Cp(p, T, patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tCp.ref() +=
            phasei().boundaryField()[patchi]*phasei().thermo().Cp(p, T, patchi);
    }

    return tCp;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::Cv() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> tCv(phasei()*phasei().thermo().Cv());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tCv.ref() += phasei()*phasei().thermo().Cv();
    }

    return tCv;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::Cv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> tCv
    (
        phasei().boundaryField()[patchi]*phasei().thermo().Cv(p, T, patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tCv.ref() +=
            phasei().boundaryField()[patchi]*phasei().thermo().Cv(p, T, patchi);
    }

    return tCv;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::gamma() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> tgamma(phasei()*phasei().thermo().gamma());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tgamma.ref() += phasei()*phasei().thermo().gamma();
    }

    return tgamma;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::gamma
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> tgamma
    (
        phasei().boundaryField()[patchi]*phasei().thermo().gamma(p, T, patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tgamma.ref() +=
            phasei().boundaryField()[patchi]
           *phasei().thermo().gamma(p, T, patchi);
    }

    return tgamma;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::Cpv() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> tCpv(phasei()*phasei().thermo().Cpv());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tCpv.ref() += phasei()*phasei().thermo().Cpv();
    }

    return tCpv;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::Cpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> tCpv
    (
        phasei().boundaryField()[patchi]*phasei().thermo().Cpv(p, T, patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tCpv.ref() +=
            phasei().boundaryField()[patchi]
           *phasei().thermo().Cpv(p, T, patchi);
    }

    return tCpv;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::CpByCpv() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> tCpByCpv(phasei()*phasei().thermo().CpByCpv());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tCpByCpv.ref() += phasei()*phasei().thermo().CpByCpv();
    }

    return tCpByCpv;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::CpByCpv
(
    const scalarField& p,
    const scalarField& T,
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> tCpByCpv
    (
        phasei().boundaryField()[patchi]*phasei().thermo().CpByCpv(p, T, patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tCpByCpv.ref() +=
            phasei().boundaryField()[patchi]
           *phasei().thermo().CpByCpv(p, T, patchi);
    }

    return tCpByCpv;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::W() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> tW(phasei()*phasei().thermo().W());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tW.ref() += phasei()*phasei().thermo().W();
    }

    return tW;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::nu() const
{
    return mu()/rho();
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::nu
(
    const label patchi
) const
{
    return mu(patchi)/rho(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::kappa() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> tkappa(phasei()*phasei().thermo().kappa());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tkappa.ref() += phasei()*phasei().thermo().kappa();
    }

    return tkappa;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::kappa
(
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> tkappa
    (
        phasei().boundaryField()[patchi]*phasei().thermo().kappa(patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tkappa.ref() +=
            phasei().boundaryField()[patchi]*phasei().thermo().kappa(patchi);
    }

    return tkappa;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::alphahe() const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<volScalarField> talphaEff(phasei()*phasei().thermo().alphahe());

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        talphaEff.ref() += phasei()*phasei().thermo().alphahe();
    }

    return talphaEff;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::alphahe
(
    const label patchi
) const
{
    PtrDictionary<phaseModel>::const_iterator phasei = phases_.begin();

    tmp<scalarField> talphaEff
    (
        phasei().boundaryField()[patchi]
       *phasei().thermo().alphahe(patchi)
    );

    for (++phasei; phasei != phases_.end(); ++phasei)
    {
        talphaEff.ref() +=
            phasei().boundaryField()[patchi]
           *phasei().thermo().alphahe(patchi);
    }

    return talphaEff;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::kappaEff
(
    const volScalarField& alphat
) const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> tkappaEff(phasei()*phasei().thermo().kappaEff(alphat));

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tkappaEff.ref() += phasei()*phasei().thermo().kappaEff(alphat);
    }

    return tkappaEff;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> tkappaEff
    (
        phasei().boundaryField()[patchi]
       *phasei().thermo().kappaEff(alphat, patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        tkappaEff.ref() +=
            phasei().boundaryField()[patchi]
           *phasei().thermo().kappaEff(alphat, patchi);
    }

    return tkappaEff;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::alphaEff
(
    const volScalarField& alphat
) const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> talphaEff(phasei()*phasei().thermo().alphaEff(alphat));

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        talphaEff.ref() += phasei()*phasei().thermo().alphaEff(alphat);
    }

    return talphaEff;
}


Foam::tmp<Foam::scalarField> Foam::multiphaseMixtureThermo::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    auto phasei = phases_.cbegin();

    tmp<scalarField> talphaEff
    (
        phasei().boundaryField()[patchi]
       *phasei().thermo().alphaEff(alphat, patchi)
    );

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        talphaEff.ref() +=
            phasei().boundaryField()[patchi]
           *phasei().thermo().alphaEff(alphat, patchi);
    }

    return talphaEff;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::rCv() const
{
    auto phasei = phases_.cbegin();

    tmp<volScalarField> trCv(phasei()/phasei().thermo().Cv());

    for (++phasei; phasei != phases_.cend(); ++phasei)
    {
        trCv.ref() += phasei()/phasei().thermo().Cv();
    }

    return trCv;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiphaseMixtureThermo::surfaceTensionForce() const
{
    tmp<surfaceScalarField> tstf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "surfaceTensionForce",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(dimensionSet(1, -2, -2, 0, 0), Zero)
        )
    );

    surfaceScalarField& stf = tstf.ref();
    stf.setOriented();

    forAllConstIters(phases_, phase1)
    {
        const phaseModel& alpha1 = *phase1;

        auto phase2 = phase1;

        for (++phase2; phase2 != phases_.cend(); ++phase2)
        {
            const phaseModel& alpha2 = *phase2;

            auto sigma = sigmas_.cfind(interfacePair(alpha1, alpha2));

            if (!sigma.good())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << " in list of sigma values"
                    << exit(FatalError);
            }

            stf += dimensionedScalar("sigma", dimSigma_, *sigma)
               *fvc::interpolate(K(alpha1, alpha2))*
                (
                    fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                  - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                );
        }
    }

    return tstf;
}


Foam::tmp<Foam::volVectorField>
Foam::multiphaseMixtureThermo::MarangoniForce(
const volScalarField& Temperature
) const
{
    tmp<volVectorField> tsMf
    (
        new volVectorField
        (
            IOobject
            (
                "MarangoniForce",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector
            (
                "MarangoniForce",
                dimensionSet(1, -2, -2, 0, 0),
                vector::zero
            )
        )
    );

    volVectorField& sMf = tsMf.ref();
    // sMf.setOriented();

    forAllConstIters(phases_, iter1)
    {
        const phaseModel& alpha1 = iter1();

        auto iter2 = iter1;

        for (++iter2; iter2 != phases_.cend(); ++iter2)
        {
            const phaseModel& alpha2 = iter2();

            auto dsigmadT = dsigmadT_.cfind(interfacePair(alpha1, alpha2));

            if (!dsigmadT.found())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << " in list of dsigmadT values"
                    << exit(FatalError);
            }

            // dsigmadTTable::const_iterator dsigmadT =
            //     dsigmadT_.find(interfacePair(alpha1, alpha2));

            // if (dsigmadT == dsigmadT_.end())
            // {
            //     FatalErrorInFunction
            //         << "Cannot find interface " << interfacePair(alpha1, alpha2)
            //         << " in list of dsigmadT values"
            //         << exit(FatalError);
            // }

                // Cell gradient of alpha
        const volVectorField gradAlpha
        (
            alpha2*fvc::grad(alpha1) - alpha1*fvc::grad(alpha2)
        );

        const volVectorField nHatM(gradAlpha/(mag(gradAlpha) + deltaN_));

        const volVectorField gradT(fvc::grad(Temperature));



        sMf -=
            dimensionedScalar("dsigmadT", dimdSigmadT_, dsigmadT())
           *(gradT-(nHatM*(nHatM & gradT)))
           *mag(gradAlpha);
        }
    }

    return tsMf;
}



Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::solve
(
    volScalarField *massdotterm
)
{
    tmp<volScalarField> tPCR
    (
        new volScalarField
        (
            IOobject
            (
                "PhaseChangeRate",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar
            (
                "PhaseChangeRate",
                dimensionSet(0, 0, -1, 0, 0),
                0.0
            )
        )
    );

    volScalarField& PCR = tPCR.ref();

    *massdotterm*=0.0;



    const Time& runTime = mesh_.time();

    const dictionary& alphaControls = mesh_.solverDict("alpha");
    label nAlphaSubCycles(alphaControls.get<label>("nAlphaSubCycles"));
    // scalar cAlpha(alphaControls.get<scalar>("cAlpha"));

    volScalarField& alpha = phases_.first();

    if (nAlphaSubCycles > 1)
    {
        surfaceScalarField rhoPhiSum(0.0*rhoPhi_);
        dimensionedScalar totalDeltaT = runTime.deltaT();

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            PCR = solveAlphas(massdotterm);
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi_;
        }

        rhoPhi_ = rhoPhiSum;
    }
    else
    {
        PCR = solveAlphas(massdotterm);
        // solveAlphas(cAlpha);
    }

    return tPCR;
}


Foam::tmp<Foam::surfaceVectorField> Foam::multiphaseMixtureThermo::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    /*
    // Cell gradient of alpha
    volVectorField gradAlpha =
        alpha2*fvc::grad(alpha1) - alpha1*fvc::grad(alpha2);

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf = fvc::interpolate(gradAlpha);
    */

    surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alpha2))
    );

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN_);
}


Foam::tmp<Foam::surfaceScalarField> Foam::multiphaseMixtureThermo::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    // Face unit interface normal flux
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}


// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::multiphaseMixtureThermo::correctContactAngle
(
    const phaseModel& alpha1,
    const phaseModel& alpha2,
    surfaceVectorField::Boundary& nHatb
) const
{
    const volScalarField::Boundary& gbf
        = alpha1.boundaryField();

    const fvBoundaryMesh& boundary = mesh_.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(gbf[patchi]))
        {
            const alphaContactAngleFvPatchScalarField& acap =
                refCast<const alphaContactAngleFvPatchScalarField>(gbf[patchi]);

            vectorField& nHatPatch = nHatb[patchi];

            vectorField AfHatPatch
            (
                mesh_.Sf().boundaryField()[patchi]
               /mesh_.magSf().boundaryField()[patchi]
            );

            const auto tp =
                acap.thetaProps().cfind(interfacePair(alpha1, alpha2));

            if (!tp.good())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << "\n    in table of theta properties for patch "
                    << acap.patch().name()
                    << exit(FatalError);
            }

            const bool matched = (tp.key().first() == alpha1.name());

            const scalar theta0 = degToRad(tp().theta0(matched));
            scalarField theta(boundary[patchi].size(), theta0);

            const scalar uTheta = tp().uTheta();

            // Calculate the dynamic contact angle if required
            if (uTheta > SMALL)
            {
                const scalar thetaA = degToRad(tp().thetaA(matched));
                const scalar thetaR = degToRad(tp().thetaR(matched));

                // Calculated the component of the velocity parallel to the wall
                vectorField Uwall
                (
                    U_.boundaryField()[patchi].patchInternalField()
                  - U_.boundaryField()[patchi]
                );
                Uwall -= (AfHatPatch & Uwall)*AfHatPatch;

                // Find the direction of the interface parallel to the wall
                vectorField nWall
                (
                    nHatPatch - (AfHatPatch & nHatPatch)*AfHatPatch
                );

                // Normalise nWall
                nWall /= (mag(nWall) + SMALL);

                // Calculate Uwall resolved normal to the interface parallel to
                // the interface
                scalarField uwall(nWall & Uwall);

                theta += (thetaA - thetaR)*tanh(uwall/uTheta);
            }


            // Reset nHatPatch to correspond to the contact angle

            scalarField a12(nHatPatch & AfHatPatch);

            scalarField b1(cos(theta));

            scalarField b2(nHatPatch.size());

            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatPatch = a*AfHatPatch + b*nHatPatch;

            nHatPatch /= (mag(nHatPatch) + deltaN_.value());
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::K
(
    const phaseModel& alpha1,
    const phaseModel& alpha2
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(alpha1, alpha2);

    correctContactAngle(alpha1, alpha2, tnHatfv.ref().boundaryFieldRef());

    // Simple expression for curvature
    return -fvc::div(tnHatfv & mesh_.Sf());
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseMixtureThermo::nearInterface() const
{
    auto tnearInt = volScalarField::New
    (
        "nearInterface",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    for (const phaseModel& phase : phases_)
    {
        tnearInt.ref() =
            max(tnearInt(), pos0(phase - 0.01)*pos0(0.99 - phase));
    }

    return tnearInt;
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::solveAlphas
(
    volScalarField *massdotterm
)
{
    //static label nSolves(-1);
    //++nSolves;

    const word alphaScheme("div(phi,alpha)");
    const word alpharScheme("div(phirb,alpha)");

    // surfaceScalarField phic(mag(phi_/mesh_.magSf()));
    // phic = min(cAlpha*phic, max(phic));


    tmp<volScalarField> tPCR
    (
        new volScalarField
        (
            IOobject
            (
                "PhaseChangeRate",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar
            (
                "PhaseChangeRate",
                dimensionSet(0, 0, -1, 0, 0),
                0.0
            )
        )
    );

    volScalarField& PCR = tPCR.ref();

    *massdotterm *= 0.0;

    PtrList<volScalarField> Sps(phases_.size());
    PtrList<volScalarField> Sus(phases_.size());

    volScalarField condensate
    (
        IOobject
        (
            "condensate",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("condensate", dimensionSet(0, 0, 0, 0, 0), 1.0)
    );

    volScalarField alphasum
    (
        IOobject
        (
            "alphasum",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("alphasum", dimensionSet(0, 0, 0, 0, 0), 0.0)
    );






    IOdictionary phasedictionary
    (
        IOobject
        (
        "thermophysicalProperties",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ_IF_MODIFIED
        )
    );




    //int fluiditer = 0;
    for (phaseModel& alpha : phases_)
    {
        
        if (alpha.isGaseous())
        {
            condensate -= alpha;
        }

        alphasum += alpha;
        //fluiditer++;
    }








    PtrList<surfaceScalarField> alphaPhiCorrs(phases_.size());

    int phasei = 0;
    for (phaseModel& alpha : phases_)
    {
        alphaPhiCorrs.set
        (
            phasei,
            new surfaceScalarField
            (
                phi_.name() + alpha.name(),
                fvc::flux
                (
                    phi_,
                    alpha,
                    alphaScheme
                )
            )
        );





                volScalarField coefffield
    (
        IOobject
        (
            "coefffield",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("coefffield", dimensionSet(0, 2, -1, 0, 0), 0)
    );



        volScalarField evaprate
    (
        IOobject
        (
            "evaprate",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("evaprate", dimensionSet(0, 0, -1, 0, 0), 0)
    );

    volScalarField condrate
    (
        IOobject
        (
            "condrate",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("condrate", dimensionSet(0, 0, -1, 0, 0), 0)
    );





            volScalarField alphagen
    (
        IOobject
        (
            "alphagen",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("alphagen", dimensionSet(0, 0, -1, 0, 0), 0)
    );

    volScalarField divU
    (
        IOobject
        (
            "divU",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("divU", dimensionSet(0, 0, -1, 0, 0), 0)
    );


    volScalarField rhogen
    (
        IOobject
        (
            "rhogen",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("rhogen", dimensionSet(1, -3, 0, 0, 0), 0)
    );

        volScalarField massgen
    (
        IOobject
        (
            "massgen",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        //dimensionedScalar("massgen", dimensionSet(1, -1, -3, 0, 0), 0)
         dimensionedScalar("massgen", dimensionSet(1, -3, -1, 1, 0), 0)
    );




        surfaceScalarField& alphaPhiCorr = alphaPhiCorrs[phasei];

        for (phaseModel& alpha2 : phases_)
        {
            if (&alpha2 == &alpha) continue;

            // surfaceScalarField phir(phic*nHatf(alpha, alpha2));



            surfaceScalarField phic
                (
                    (mag(phi_))/mesh_.magSf()
                );

                surfaceScalarField phir(0.0*phic*nHatf(condensate, (1.0-condensate)));//JUST TO INITIALISE TO ZERO - must be a neater way

                surfaceScalarField phirD(0.0*mag(phi_/mesh_.magSf())*nHatf(condensate, (1.0-condensate)));//setto zero
                phirD*=0.0;

                scalarCoeffSymmCTable::const_iterator cAlpha
            (
                cAlphas_.find(interfacePair(alpha, alpha2))
            );


            if (cAlpha != cAlphas_.end())
            {


                surfaceScalarField phic2
                (
                    (mag(phi_) + mag(phir))/mesh_.magSf()
                );



                phir +=(min(cAlpha()*phic2, max(phic2))*(nHatf(alpha, alpha2)));


            }

            scalarCoeffSymmDTable::const_iterator dAlpha
            (
                dAlphas_.find(interfacePair(alpha, alpha2))
            );

            if (dAlpha != dAlphas_.end())
            {
            dimensionedScalar valdiff("valdiff",dimdiff_,dAlpha());

            coefffield=(/*epsilon1*/valdiff);



            phirD-= fvc::interpolate(coefffield)*mesh_.magSf()*((fvc::interpolate(alpha2)*fvc::snGrad(alpha))-(fvc::interpolate(alpha)*fvc::snGrad(alpha2)));;



            }



        dimensionedScalar maxrate( "maxrate", dimensionSet(0,0,-1,0,0,0,0), 0.5/mesh_.time().deltaT().value() );//something like rate
        dimensionedScalar gasconstant("gasconstant",dimensionSet(1, 2, -2, -1, -1),scalar(8.314));//1 2 -2 -1 -1 proper units


        dimensionedScalar P0("P0",dimensionSet(1, -1, -2, 0, 0),phasedictionary);//1 2 -2 -1 -1 proper units

        // Info<<alpha.name()<<"\t mol weight \t"<<alpha.thermo().W();

        boilTable::const_iterator boilT
            (
                boils_.find(interfacePair(alpha, alpha2))
            );

            if(boilT != boils_.end()){// return boiling temperature for pair

            dimensionedScalar pair_boil_T("pair_boil_T",dimBoil_,boilT());


            ////// iterator to look up latent heat gas of phase paire
            LatentHeatGasTable::const_iterator LHG
            (
                LatentHeatGass_.find(interfacePair(alpha, alpha2))
            );
            dimensionedScalar pair_LHG("pair_LHG",dimLatentHeatGas_,LHG());

            // Info<<"TEST_LHG"<<alpha.name()<<"\t"<<alpha2.name()<<"\t"<<pair_LHG<<endl;
////// iterator to look up latent heat gas of phase paire



    volScalarField Psat
    (
        IOobject
        (
            "Psat",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("Psat", dimensionSet(1, -1, -2, 0, 0), 0.0)
    );

    Psat = P0*Foam::exp((((alpha.thermo().W()/1000.0)*pair_LHG)/(pair_boil_T*gasconstant))*(1.0-(pair_boil_T/T_)));

    // Info<<Psat<<endl;


    volScalarField evapcoefffield
    (
        IOobject
        (
            "evapcoefffield",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("evapcoefffield", dimensionSet(0, 0, 0, 0, 0), 0.0),
            zeroGradientFvPatchScalarField::typeName
    );



// if(p_>Psat){
// evapcoefffield=0.0;//condensation
// }
// else{
// evapcoefffield=1.0;//evaporation
// }


   forAll( mesh_.C(), celli)
{
if(p_[celli]>Psat[celli]){
evapcoefffield[celli]=0.0;
}
else{
evapcoefffield[celli]=1.0;//(1.0/(alpha2[celli]));
// Info<<evapcoefffield[celli]<<endl;
}
}

// Info<<evapcoefffield<<endl;

dimensionedScalar int_thickness("interface_thickness",dimensionSet(0, 1, 0, 0, 0),phasedictionary);


if(alpha2.name()==(alpha.name()+"vapour")){//alpha2 is vapour , alpha1 is liquid
volScalarField liqdensity(alpha.thermo().rho());
volScalarField vapdensity(alpha2.thermo().rho());
Info<<"Liquid-Vapour State Transition: (Liquid,Vapour): ("<<alpha.name()<<","<<alpha2.name()<<")  Boiling Temperature: "<<pair_boil_T<<"K"<<endl;

                        // Info<<"HERE1"<<endl;
                        //NEW Compressible cond and evap rates
                        evaprate = 1.0*Foam::sqrt((alpha.thermo().W()/1000.0)/(2.0*M_PI*gasconstant*T_))*(1.0/(int_thickness*liqdensity))*(Psat-p_);//*****************;
                        condrate = -max(condensate,1e-6)*Foam::sqrt((alpha2.thermo().W()/1000.0)/(2.0*M_PI*gasconstant*T_))*(1.0/(int_thickness*vapdensity))*(Psat-p_);//**********************;
                        //NEW Compressible cond and evap rates
                        // Info<<"HERE2"<<endl;
                        // evaprate=ratedummy;//L_Lee*(alpha.thermo().rho()/alpha2.thermo().rho())*((T_-pair_boil_T)/pair_boil_T);//Lee

                        // condrate=ratedummy;//L_Lee*(alpha2.thermo().rho()/alpha.thermo().rho())*((pair_boil_T-T_)/pair_boil_T);//Lee


                    // alphagen=min(condrate,maxrate)*(1.0-evapcoefffield)*alpha2*(alpha2.thermo().rho()/alpha.thermo().rho());//0.1;//0.1;//*(alpha2.thermo().rho()/alpha.thermo().rho());
                    alphagen=((min(condrate,maxrate)*(1.0-evapcoefffield)*alpha2*(alpha2.thermo().rho()/alpha.thermo().rho()))-(min(evaprate,maxrate)*evapcoefffield*alpha*(alpha2.thermo().rho()/alpha.thermo().rho())));
                        // divU = -condrate*alpha2*10.0;
                         massgen=pair_LHG*(((alpha.thermo().rho()/alpha.thermo().Cv())*min(condrate,maxrate)*(1.0-evapcoefffield)*alpha2*(alpha2.thermo().rho()/alpha.thermo().rho()))-((alpha2.thermo().rho()/alpha2.thermo().Cv())*min(evaprate,maxrate)*evapcoefffield*alpha*(alpha2.thermo().rho()/alpha.thermo().rho())));


            }

            else if(alpha.name()==(alpha2.name()+"vapour")){//alpha is vapour to phase2
            volScalarField liqdensity(alpha2.thermo().rho());
            volScalarField vapdensity(alpha.thermo().rho());

                    //NEW Compressible cond and evap rates
                    // Info<<"HERE3"<<endl;
                        evaprate = 1.0*Foam::sqrt((alpha2.thermo().W()/1000.0)/(2.0*M_PI*gasconstant*T_))*(1.0/(int_thickness*liqdensity))*(Psat-p_);//*****************;
                        condrate = -max(condensate,1e-6)*Foam::sqrt((alpha.thermo().W()/1000.0)/(2.0*M_PI*gasconstant*T_))*(1.0/(int_thickness*vapdensity))*(Psat-p_);//**********************;
                        //NEW Compressible cond and evap rates
// Info<<"HERE4"<<endl;




                        // evaprate=ratedummy;//L_Lee*(alpha2.thermo().rho()/alpha.thermo().rho())*((T_-pair_boil_T)/pair_boil_T);

                        // condrate=ratedummy;//L_Lee*(alpha.thermo().rho()/alpha2.thermo().rho())*((pair_boil_T-T_)/pair_boil_T);


                        // alphagen=((min(evaprate,maxrate)*evapcoefffield*alpha2)-(min(condrate,maxrate)*(1.0-evapcoefffield)*alpha));
                        alphagen=((min(evaprate,maxrate)*evapcoefffield*alpha2)-(min(condrate,maxrate)*(1.0-evapcoefffield)*alpha));//
                    //  divU = -condrate*alpha*1.0;//alphagen;
                        // evapcondrho= -(fvc::interpolate(alpha)*100.0/*+fvc::interpolate(alpha2)*10.0*/)*rhoone;
                        massgen=pair_LHG*(((alpha.thermo().rho()/alpha.thermo().Cv())*min(evaprate,maxrate)*evapcoefffield*alpha2)-((alpha2.thermo().rho()/alpha2.thermo().Cv())*min(condrate,maxrate)*(1.0-evapcoefffield)*alpha));



            }
            }


// Info<<"HERE4"<<endl;




            alphaPhiCorr += fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha,
                alpharScheme
            )
            +phirD
            ;
        }


// Info<<"HERE5"<<endl;

    PCR+=alphagen;//alphagen;//divU;//

    // PCR-=(1.0-(alphasum))*(10e-2/mesh_.time().deltaT());

    *massdotterm+=massgen;


        Sps.set(phasei,new volScalarField(divU));//
        volScalarField& Sp=Sps[phasei];
        Sus.set(phasei,new volScalarField(alphagen));//
        volScalarField& Su=Sus[phasei];

// Info<<"HERE6"<<endl;


        MULES::limit
        (
            1.0/mesh_.time().deltaT().value(),
            geometricOneField(),
            alpha,
            phi_,
            alphaPhiCorr,
            Sp,//zeroField(),
            Su,//alphalap, //zeroField()
            oneField(),
            zeroField(),
            true//NEEDS TO BE TRUE
        );

        ++phasei;
    }

    MULES::limitSum(alphaPhiCorrs);

    rhoPhi_ = dimensionedScalar(dimensionSet(1, 0, -1, 0, 0), Zero);

    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );


    volScalarField divU(fvc::div(fvc::absolute(phi_, U_)));


    phasei = 0;

    for (phaseModel& alpha : phases_)
    {
        surfaceScalarField& alphaPhi = alphaPhiCorrs[phasei];
        alphaPhi += upwind<scalar>(mesh_, phi_).flux(alpha);

        volScalarField::Internal Sp
        (
            IOobject
            (
                "Sp",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(alpha.dgdt().dimensions(), Zero)
        );

        volScalarField::Internal Su
        (
            IOobject
            (
                "Su",
                mesh_.time().timeName(),
                mesh_
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            divU*min(alpha, scalar(1))
            +
            Sus[phasei]
            -PCR*(min(alpha, scalar(1)))//double counting divU - need to seperate out div(U) due to phase change
        );

        {
            const scalarField& dgdt = alpha.dgdt();

            forAll(dgdt, celli)
            {
                if (dgdt[celli] < 0.0 && alpha[celli] > 0.0)
                {
                    Sp[celli] += dgdt[celli]*alpha[celli];
                    Su[celli] -= dgdt[celli]*alpha[celli];
                }
                else if (dgdt[celli] > 0.0 && alpha[celli] < 1.0)
                {
                    Sp[celli] -= dgdt[celli]*(1.0 - alpha[celli]);
                }
            }
        }

        for (const phaseModel& alpha2 : phases_)
        {
            if (&alpha2 == &alpha) continue;

            const scalarField& dgdt2 = alpha2.dgdt();

            forAll(dgdt2, celli)
            {
                if (dgdt2[celli] > 0.0 && alpha2[celli] < 1.0)
                {
                    Sp[celli] -= dgdt2[celli]*(1.0 - alpha2[celli]);
                    Su[celli] += dgdt2[celli]*alpha[celli];
                }
                else if (dgdt2[celli] < 0.0 && alpha2[celli] > 0.0)
                {
                    Sp[celli] += dgdt2[celli]*alpha2[celli];
                }
            }
        }

        MULES::explicitSolve
        (
            geometricOneField(),
            alpha,
            alphaPhi,
            Sp,
            Su
        );

        rhoPhi_ += fvc::interpolate(alpha.thermo().rho())*alphaPhi;

        Info<< alpha.name() << " volume fraction, min, max = "
            << alpha.weightedAverage(mesh_.V()).value()
            << ' ' << min(alpha).value()
            << ' ' << max(alpha).value()
            << endl;

        sumAlpha += alpha;

        ++phasei;
    }

    Info<< "Phase-sum volume fraction, min, max = "
        << sumAlpha.weightedAverage(mesh_.V()).value()
        << ' ' << min(sumAlpha).value()
        << ' ' << max(sumAlpha).value()
        << endl;










    // Calculate and print total mass
    Info << "Phase pair (liquid + vapor) masses in domain:" << endl;
    
    HashTable<bool> processedPhases;
    
    for (const phaseModel& phase : phases_)
    {
        if (processedPhases.found(phase.name())) continue;
        
        word vaporName = phase.name() + "vapour";
        bool hasVaporPair = false;
        scalar liquidMass = 0.0;
        scalar vaporMass = 0.0;
        
        const volScalarField phaseMassLiq(phase*phase.thermo().rho());
        liquidMass = gSum(phaseMassLiq.primitiveField()*mesh_.V().field());
        
        for (const phaseModel& phase2 : phases_)
        {
            if (phase2.name() == vaporName)
            {
                hasVaporPair = true;
                const volScalarField phaseMassVap(phase2*phase2.thermo().rho());
                vaporMass = gSum(phaseMassVap.primitiveField()*mesh_.V().field());
                processedPhases.insert(phase2.name(), true);
                break;
            }
        }
        
        if (hasVaporPair)
        {
            scalar totalPairMass = liquidMass + vaporMass;
            Info << "    " << phase.name() << " + " << vaporName << ": " 
                 << totalPairMass << " kg (liquid: " << liquidMass 
                 << " kg, vapor: " << vaporMass << " kg)" << endl;
            processedPhases.insert(phase.name(), true);
        }
        else if (phase.name() != "air" && phase.name().find("vapour") == string::npos)
        {
            Info << "    " << phase.name() << " (no vapor): " << liquidMass << " kg" << endl;
            processedPhases.insert(phase.name(), true);
        }
        else if (phase.name() == "air")
        {
            Info << "    " << phase.name() << ": " << liquidMass << " kg" << endl;
            processedPhases.insert(phase.name(), true);
        }
    }










    calcAlphas();

    return tPCR;
}


// ************************************************************************* //
