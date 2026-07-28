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
#include "PstreamReduceOps.H"   // reduce()/sumOp for the clamp counters
#include "subCycle.H"
#include "MULES.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDdt.H"
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
    dAlphas_(lookup("interfaceDiffusion")),
    dimdiff_(0, 2, -1, 0, 0),
    deltaN_
    (
        "deltaN",
        1e-8/cbrt(average(mesh_.V()))
    ),

    // ----- conserved-partial-mass scheme state (ADDED) ---------------------
    //  Declaration order in the header is: ... deltaN_, c_, c0_, beta0_,
    //  curTimeIndex_  -- this initialiser list must match that order.
    c_(),       // sized + seeded in the body (needs the densities from correct())
    c0_(),      // sized + seeded in the body
    beta0_
    (
        IOobject
        (
            "beta0",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    ),
    curTimeIndex_(-1)
{
    rhoPhi_.setOriented();
    calcAlphas();
    alphas_.write();
    correct();                  // therm    o (phase densities) are valid after this

         c_.setSize(phases_.size());
    c0_.setSize(phases_.size());

    label phasei = 0;
    for (phaseModel& phase : phases_)
    {
        {
            IOobject cIO
            (
                "c." + phase.name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            );

            if (cIO.typeHeaderOk<volScalarField>(true))
            {
                c_.set
                (
                    phasei,
                    new volScalarField
                    (
                        IOobject
                        (
                            "c." + phase.name(),
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::MUST_READ,
                            IOobject::AUTO_WRITE
                        ),
                        mesh_
                    )
                );
            }
            else
            {
                // Fresh start: derive once from alpha*rho (EXACt here).
                c_.set
                (
                    phasei,
                    new volScalarField(cIO, phase*phase.thermo().rho())
                );
            }
        }

        c0_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    "c0." + phase.name(),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                c_[phasei]
            )
        );

        ++phasei;
    }

    beta0_ == dimensionedScalar("zero", dimless, 0.0);
    for (const phaseModel& phase : phases_)
    {
        if (!phase.isGaseous())
        {
            beta0_ += phase;
        }
    }
}




// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::multiphaseMixtureThermo::restoreOldTimeValues()
{
    // Any term used or involbed  in `fvc::ddt` needs to have their `.oldTime()`
    // value reread from disk if solver restarted!

    
    rhoPhi_.oldTime();
    for (phaseModel& alpha : phases_)
    {
        volScalarField& rhoK = alpha.thermo().rho();
        rhoK.oldTime();
    }

    
    if (!(mesh_.time().timeIndex() > 1))
    {
        
        return;
    }

    label errorCount = 0;
    // Restore values if solver restart detected
    //- Restore rhoPhi oldTime
    {
        surfaceScalarField rhoPhi_oldTime
        (
            IOobject
            (
                IOobject::groupName("rhoPhi", "oldTime"),
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            rhoPhi_
        );

        if (rhoPhi_oldTime.headerOk())
        {
            rhoPhi_.oldTime() = rhoPhi_oldTime;
        }
        else
        {
            Info<< "Mixture: WARNING could not find field `"
                << "rhoPhi.oldTime`. "
                << "This will cause issues in fvc::ddt terms!"
                << "Please check your latest time directory (folder)!"
                << endl;

            errorCount++;
        }
    }


    for (phaseModel& alpha : phases_)
    {
        volScalarField rhoKoldTime
        (
            IOobject
            (
                "rho." + alpha.name() + ".oldTime",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            alpha.thermo().rho()
        );

        if (rhoKoldTime.headerOk())
        {
            volScalarField& rhoK = alpha.thermo().rho();
            rhoK.oldTime() = rhoKoldTime;
        }
        else
        {
            Info<< "Mixture: WARNING could not find field `"
                << "rho." << alpha.name() << ".oldTime`. "
                << "This will cause issues in fvc::ddt terms!"
                << "Please check your latest time directory (folder)!"
                << endl;

            errorCount++;
        }
    }

    if (errorCount)
    {
        Info<< "Mixture WARNING: " << errorCount << " fields' oldTime values were not "
            << "restored! This will cause issues within any fvc::ddt terms!"
            << "RECOMMENDED: please inspect your latest time folder!"
            << endl;
    }
    else
    {
        Info<< "Mixture: All relevant fields' oldTime values were restored using "
            << "file reads. All good!"
            << endl;
    }
}


void Foam::multiphaseMixtureThermo::writeOldTimeValues()
{
    if (mesh_.time().writeTime())
    {
        surfaceScalarField rhoPhi_oldTime
        (
            IOobject
            (
                IOobject::groupName("rhoPhi", "oldTime"),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            rhoPhi_.oldTime()
        );

        for (phaseModel& alpha : phases_)
        {
            volScalarField& rhoK = alpha.thermo().rho();

            volScalarField rhoKoldTime
            (
                IOobject
                (
                    "rho." + alpha.name() + ".oldTime",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                rhoK.oldTime()
            );

            rhoKoldTime.write();
        }
    }
}


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


// void Foam::multiphaseMixtureThermo::correctRho(const volScalarField& dp)
// {
//     for (phaseModel& phase : phases_)
//     {
//         phase.thermo().rho() += phase.thermo().psi()*dp;
//     }
// }

void Foam::multiphaseMixtureThermo::correctRho(const volScalarField& dp)
{
    for (phaseModel& phase : phases_)
    {
        volScalarField& rhoPh = phase.thermo().rho();
        rhoPh += phase.thermo().psi()*dp;
        rhoPh = max(rhoPh, dimensionedScalar("rhoMinEOS", dimDensity, 1e-4));
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
    volScalarField* massdotterm,
    volScalarField* vDotPptr
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


    const dictionary& alphaControls = mesh_.solverDict("alpha");
    const label nAlphaSubCycles
    (
        alphaControls.getOrDefault<label>("nAlphaSubCycles", 1)
    );

    if (nAlphaSubCycles > 1)
    {
        FatalErrorInFunction
            << "nAlphaSubCycles = " << nAlphaSubCycles
            << ", but this solver's conservative species scheme"
            << " requires nAlphaSubCycles == 1 (sub -cycling breaks this"
            << " )." << nl
            << "Set nAlphaSubCycles 1 in fvSolution. We'll try and allow sub-stepping in the future but for now chill out and just set it to 1"
            << exit(FatalError);
    }

    PCR = solveAlphas(massdotterm, vDotPptr);

    writeOldTimeValues();

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

void Foam::multiphaseMixtureThermo::correctContactAngle//from olesr bersion
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


//  multiphaseMixtureThermo::solveAlphas   --   rewriting this to be fully mass conservative through a kind of heirarchical volume of fluid aproach
//                                           there is one condensed field that is limited by MULES to keep the interface sharp, but then explicit partial density transport
//
// conserved-partial-mass formulation
//   conserved variables are the per-species PARTIAL MASSES
//
//        c_i = alpha_i * rho_i        units are  [kg/m^3]
//

Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::solveAlphas
(
    volScalarField* massdotterm,
    volScalarField* vDotPptr
)
{
    const word alphaScheme("div(phi,alpha)");
    const word alpharScheme("div(phirb,alpha)");

    const Time& runTime   = mesh_.time();
    const scalar dt       = runTime.deltaT().value();
    const scalar rDeltaT  = 1.0/dt;

    // ----- return field: phase-change volumetric rate (-> vDot) -------------
    tmp<volScalarField> tPCR
    (
        new volScalarField
        (
            IOobject("PhaseChangeRate", runTime.timeName(), mesh_),
            mesh_,
            dimensionedScalar("PhaseChangeRate", dimensionSet(0, 0, -1, 0, 0), 0.0)
        )
    );
    volScalarField& PCR = tPCR.ref();

    *massdotterm *= 0.0;

    if (vDotPptr)
    {
        *vDotPptr *= 0.0;
    }


    const dimensionedScalar Tmin("Tmin", dimTemperature, 300.0);// just for the FPE that is happening, just a catch shouldnt change physics especially is timestep is below Co etc
    const volScalarField Tsafe(max(T_, Tmin)); //  NEW ADDED BY TOM July to stop blow up doesnt affect energy juts rates of transfer

    IOdictionary phasedictionary
    (
        IOobject
        (
            "thermophysicalProperties",
            runTime.constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED
        )
    );
    const dimensionedScalar gasconstant
        ("gasconstant", dimensionSet(1, 2, -2, -1, -1), scalar(8.314));
    const dimensionedScalar P0
        ("P0", dimensionSet(1, -1, -2, 0, 0), phasedictionary);
    const dimensionedScalar int_thickness
        ("interface_thickness", dimensionSet(0, 1, 0, 0, 0), phasedictionary);

    // Per-step volumetric dilatation limit for the phase-change
    const scalar etaVol
    (
        phasedictionary.getOrDefault<scalar>("phaseChangeVolLimit", 0.1)
    );

    // Hertz-Knudsen evaporation/condensation accomodation coefficient, //  NEW ADDED BY TOM July
    // applied symmetrically to both directions and to the implicit
    // pressure-coupling coefficient kP.  ~0.8-1.0 for metals.
    // Roman, this is a new part where we basically partition the rates rather than using a pressure per material
    const scalar accommodationCoeff
    (
        phasedictionary.getOrDefault<scalar>("accommodationCoeff", 1.0)
    );

    const label nPhases = phases_.size();


    //  Phase bookeeping: index map + condensed/gas fields.
    //  phaseModel is-a volScalarField; isGaseous() distinguishes vapour/air.
     //  NEW ADDED BY TOM July

    HashTable<label> name2idx(2*nPhases);
    List<bool>  isGas(nPhases, false);
    {
        label i = 0;
        for (const phaseModel& ph : phases_)
        {
            name2idx.insert(ph.name(), i);
            isGas[i] = ph.isGaseous();
            ++i;
        }
    }


    volScalarField condensate
    (
        IOobject("condensate", runTime.timeName(), mesh_),
        mesh_,
        dimensionedScalar("condensate", dimless, 1.0)
    );
    for (const phaseModel& ph : phases_)
    {
        if (ph.isGaseous()) condensate -= ph;
    }
    // numerically pin to [0,1] for the thing gooing through the change value (MULES re-bounds)
    condensate = max(min(condensate, scalar(1)), scalar(0));

    // --- Interface gate for the HK kinetics ---------------------------------
    // r is an interface flux spread over int_thickness; ungated it fires in
    // PURE cells too, with the driving degenerating at both ends:
    //   * beta -> 0 (gas-only cavity/plume): xLiq comes from the trace
    //     condensed residue while y*p is real, so the driving is generically
    //     negative: sustained condensation-to-fog with NO liquid present.
    //     The mass clamp (0.5*c_vap/dt) only checks VAPOUR availability, so
    //     the volume sink runs continuously, and kP (nonzero there) makes the
    //     implicit pEqn term relax cavity pressure toward the fog equilibrium
    //     x*Psat/y -- typically BELOW pMin: the persistently pinned floor and
    //     the air-filled keyhole.
    //   * beta -> 1 (bulk superheated melt): yVap -> 0 so the driving becomes
    //     +x*Psat: spurious sub-surface flash boiling with no interface.
    // Gate on co-presence of both groups.  Saturates to 1 across the resolved
    // interface band (min(beta,1-beta) >= betaGate), so calibrated interface
    // rates are untouched; exactly zero in pure cells.  Supersaturated vapour
    // then simply advects until it reaches an interface (no homogeneous
    // nucleation model -- add an area-density closure if plume fog is wanted).
    const scalar betaGate
    (
        phasedictionary.getOrDefault<scalar>("phaseChangeGate", 0.01)
    );
    const volScalarField gInt
    (
        IOobject("phaseChangeGateField", runTime.timeName(), mesh_),
        min
        (
            min(condensate, scalar(1) - condensate)/betaGate,
            scalar(1)
        )
    );


    if (runTime.timeIndex() != curTimeIndex_)
    {
        curTimeIndex_ = runTime.timeIndex();
        beta0_ == condensate;
        forAll(c_, i)
        {
            c0_[i] == c_[i];
        }
    }
    forAll(c_, i)
    {
        c_[i].oldTime() == c0_[i];
    }


    //  Compressibility fields for the bounded dilatation/relative split.
    //    dgdtBar      = SUM_all   alpha_j dgdt_j     (volume-averaged)
    //    dgdtBetaSum  = SUM_cond  alpha_i dgdt_i

    volScalarField dgdtBar
    (
        IOobject("dgdtBar", runTime.timeName(), mesh_),
        mesh_,
        dimensionedScalar("dgdtBar", dimensionSet(0, 0, -1, 0, 0), 0.0)
    );
    volScalarField dgdtBetaSum
    (
        IOobject("dgdtBetaSum", runTime.timeName(), mesh_),
        mesh_,
        dimensionedScalar("dgdtBetaSum", dimensionSet(0, 0, -1, 0, 0), 0.0)
    );
    {
        label i = 0;
        for (phaseModel& ph : phases_)
        {
            dgdtBar += ph*ph.dgdt();
            if (!isGas[i]) dgdtBetaSum += ph*ph.dgdt();
            ++i;
        }
    }

    //  NEW ADDED BY TOM July
                     //  RAOULT-DALTON MULTICOMPONENT KINETICS:
    //   * x_i = liquid mole fraction of species i over the CONDENSED group
    //     
    //   * y_i = vapour mole fraction over the GASEOUS group incl. air
    //     (Dalton: the gas-side resistance is the PARTIAL pressure y_i*p,
    //     not the total p);


    //   * sigma_e = accommodation coefficient (dict `accommodationCoeff`,
    //     default 1; ~0.8-1 for clean liquid metals).


   
    const volScalarField rCvMix(this->rCv());   // mixture reciprocal Cv

    const volScalarField rhoMixField(this->rho());


    //
    //  Molar densities   n_i = max(c_i, 0)/W_i           
    //  Liquid mole fracs x_i = n_i / SUM_condensed n     
    //  Gas mole fracs    y_i = n_i / SUM_gaseous  n      

    PtrList<volScalarField> moleFrac(nPhases);
    {
        const dimensionedScalar nSmall
            ("nSmall", dimMoles/dimVolume, VSMALL);
        const dimensionedScalar cZeroML
            ("cZeroML", dimDensity, 0.0);

        volScalarField nCond
        (
            IOobject("nCond", runTime.timeName(), mesh_),
            mesh_,
            dimensionedScalar("0", dimMoles/dimVolume, 0.0)
        );
        volScalarField nGas
        (
            IOobject("nGas", runTime.timeName(), mesh_),
            mesh_,
            dimensionedScalar("0", dimMoles/dimVolume, 0.0)
        );

        PtrList<volScalarField> nMol(nPhases); //  NEW ADDED BY TOM July
        {
            label i = 0;
            for (const phaseModel& ph : phases_)
            {
                nMol.set
                (
                    i,
                    new volScalarField
                    (
                        "n." + ph.name(),
                        max(c_[i], cZeroML)/ph.thermo().W()
                    )
                );
                if (isGas[i]) nGas += nMol[i];
                else          nCond += nMol[i];
                ++i;
            }
        }

        label i = 0;
        for (const phaseModel& ph : phases_)
        {
            const volScalarField& nHost = isGas[i] ? nGas : nCond; //  NEW ADDED BY TOM July think this is pretty nice but to be checked

            moleFrac.set
            (
                i,
                new volScalarField
                (
                    (isGas[i] ? "y." : "x.") + ph.name(),
                    min
                    (
                        max(nMol[i]/max(nHost, nSmall), scalar(0)),
                        scalar(1)
                    )
                )
            );
            ++i;
        }
    }

    PtrList<volScalarField> SuMass(nPhases);     // mass source per phase [kg/m^3/s] .... beter naming than old version
    forAll(SuMass, i)
    {
        SuMass.set
        (
            i,
            new volScalarField
            (
                IOobject("SuMass" + Foam::name(i), runTime.timeName(), mesh_),
                mesh_,
                dimensionedScalar("0", dimensionSet(1, -3, -1, 0, 0), 0.0)
            )
        );
    }

    Info<< "Phase-change pairs:" << endl;
    {
        label li = 0;
        for (phaseModel& liq : phases_)
        {
            if (!isGas[li])     // liquid / condensed metal phase
            {
                const word vapName = liq.name() + "vapour";
                if (name2idx.found(vapName))
                {
                    const label vi = name2idx[vapName];

                    phaseModel& vap = phases_[vapName];

                    // pair lookups (Tboil, latent heat) from the interface tables .... gone back to this 2 table way,,,, maybe update to that new class itried where this is in one object..
                    boilTable::const_iterator boilT
                        (boils_.find(interfacePair(liq, vap)));
                    LatentHeatGasTable::const_iterator LHG
                        (LatentHeatGass_.find(interfacePair(liq, vap)));

                    if (boilT != boils_.end() && LHG != LatentHeatGass_.end())
                    {
                        const dimensionedScalar Tboil
                            ("Tboil", dimBoil_, boilT());
                        const dimensionedScalar Lvap
                            ("Lvap", dimLatentHeatGas_, LHG());

                        const volScalarField Wliq
                        (
                            IOobject("Wliq." + liq.name(), runTime.timeName(), mesh_),
                            liq.thermo().W()
                        );


                        const dimensionedScalar rhoFloorPC// dont let density go below this as it messes up the solve and crashes.. maybe try other values
                            ("rhoFloorPC", dimDensity, 1e-3);
                        const volScalarField rhoLiq(max(liq.thermo().rho(), rhoFloorPC));
                        const volScalarField rhoVap(max(vap.thermo().rho(), rhoFloorPC));

                        // Raoult-Dalton composition of this pair (ledger-based)
                        const volScalarField& xLiq = moleFrac[li];
                        const volScalarField& yVap = moleFrac[vi];

                        // Clausius-Clapeyron saturation pressure (legacy form)
                        const volScalarField Psat
                        (
                            IOobject("Psat." + liq.name(), runTime.timeName(), mesh_),
                            P0*Foam::exp
                            (
                                (((Wliq/1000.0)*Lvap)/(Tboil*gasconstant))
                               *(1.0 - (Tboil/Tsafe))
                            )
                        );

                        // unified signed mass rate  [kg/m^3/s]................. //  NEW ADDED BY TOM July
                                // Raoult-Dalton driving force: liquid-side equilibrium
                         // partial pressure x_i*Psat_i minus vapour-side partial
                        // pressure y_i*p.  
                        volScalarField r
                        (
                            IOobject("r" + liq.name(), runTime.timeName(), mesh_),
                            accommodationCoeff
                           *Foam::sqrt
                            (
                                (Wliq/1000.0)/(2.0*M_PI*gasconstant*Tsafe)
                            )
                           *(xLiq*Psat - yVap*p_)/int_thickness
                        );

                        // interface gating (see gInt construction above):
                        // applied BEFORE the clamps so unclampedI and kP see
                        // the gated kinetic rate consistently.
                        r *= gInt;

                        
                        //  Composition-consistent saturation temperature: the T//////////////// //  NEW ADDED BY TOM July
                        //  at which THIS pair's Raoult-Dalton driving force is
                        //  exactly zero,
                        //        x_i*Psat_i(Tsat) = y_i*p,/////////////// //  NEW ADDED BY TOM July PAPER REF NEEDED 
                        //  so the energy-to-saturation clamp brakes toward the
                        //  SAME fixed point as the kinetics.  Inverting
                        //        Psat = P0 exp( K (1 - Tboil/T) ),
                        //        K = (W/1000) Lvap / (Tboil R),
                        //  gives  Tsat = Tboil / ( 1 - ln( y p /(x P0) )/K ).
                        
                        const volScalarField Kfield //for above tect relation///
                        (
                            ((Wliq/1000.0)*Lvap)/(Tboil*gasconstant)
                        );
                        const volScalarField satArg
                        (
                            max
                            (
                                yVap
                               *max
                                (
                                    p_,
                                    dimensionedScalar
                                    ("pSmallSat", P0.dimensions(), 1.0)
                                )
                               /(
                                    max
                                    (
                                        xLiq,
                                        dimensionedScalar
                                        ("xSmall", dimless, 1e-12)
                                    )*P0
                                ),
                                dimensionedScalar("argSmall", dimless, 1e-15)
                            )
                        );
                        const volScalarField Tsat
                        (
                            IOobject("Tsat." + liq.name(), runTime.timeName(), mesh_),
                            Tboil
                          / max
                            (
                                1.0 - Foam::log(satArg)/Kfield,
                                dimensionedScalar("d", dimless, 0.05)
                            )
                        );


                        scalarField unclampedI(r.primitiveField().size(), 1.0);//Clamping to make sure we dont go unstavle
                        label nClampE = 0, nClampM = 0, nClampV = 0;
                        {
                          

                            scalarField& rI = r.primitiveFieldRef();
                            const scalarField& TI      = Tsafe.primitiveField();
                            const scalarField& TsatI   = Tsat.primitiveField();
                            const scalarField& rhoMixI = rhoMixField.primitiveField();
                            const scalarField& rCvMixI = rCvMix.primitiveField();
                            const scalar       LvapV   = Lvap.value();
                            const scalarField& cLiqI = c_[li].primitiveField();
                            const scalarField& cVapI = c_[vi].primitiveField();
                            const scalarField& rhoLiqI = rhoLiq.primitiveField();
                            const scalarField& rhoVapI = rhoVap.primitiveField();
                            const scalar cvSmall = 1e-30;
                            forAll(rI, celli)//LIMITING SO CODE DOESNT BLOW UP
                            {
                                scalar& rc = rI[celli];
                                const scalar rcKinetic = rc;

                                //  energy limit 
                                const scalar eFac =
                                    rhoMixI[celli]*rDeltaT
                                  / (Foam::max(rCvMixI[celli], cvSmall)*LvapV);
                                const scalar headUp =
                                    Foam::max(TsatI[celli] - TI[celli], 0.0);
                                const scalar headDn =
                                    Foam::max(TI[celli] - TsatI[celli], 0.0);
                                const scalar rMinE = -headUp*eFac;  // cond cap (<=0)
                                const scalar rMaxE =  headDn*eFac;  // evap cap (>=0)
                                if (rc > rMaxE) { rc = rMaxE; ++nClampE; }
                                if (rc < rMinE) { rc = rMinE; ++nClampE; }

                                //  mass available
                                const scalar rmaxEvap =
                                    0.5*Foam::max(cLiqI[celli], 0.0)*rDeltaT;
                                const scalar rmaxCond =
                                    0.5*Foam::max(cVapI[celli], 0.0)*rDeltaT;
                                if (rc >  rmaxEvap) { rc =  rmaxEvap; ++nClampM; } // evap
                                if (rc < -rmaxCond) { rc = -rmaxCond; ++nClampM; } // cond

                                //  VOLUMETRIC 
                                const scalar dVspec =
                                    Foam::max
                                    (
                                        1.0/rhoVapI[celli]
                                      - 1.0/rhoLiqI[celli],
                                        VSMALL
                                    );
                                const scalar rMaxV = etaVol*rDeltaT/dVspec;
                                if (rc >  rMaxV) { rc =  rMaxV; ++nClampV; }
                                if (rc < -rMaxV) { rc = -rMaxV; ++nClampV; }


                                if (rc != rcKinetic)
                                {
                                    unclampedI[celli] = 0.0;
                                }
                            }
                        }

                        // exactly-paired mass sources (same rate r)
                        SuMass[li] -= r;     // liquid loses
                        SuMass[vi] += r;     // vapour gains (paired)

                        // volumetric expansion goes to pressure equation through vDot
                        PCR += r*(1.0/rhoVap - 1.0/rhoLiq);

                        // latent-heat of vap/cond goes to TEqn through -mass_dot ..... maybe i should change the sign, but this is right !!!
                        *massdotterm += rCvMix*Lvap*r;

                        
                        //  Implicit pressure-coupling goes into vdot
                        //

                 //      vDot_pair = sigma_e*A(T)/delta
          //                 *(x_i*Psat_i - y_i*p)*dVspec,
             //     with  A(T) = sqrt((W/1000)/(2 pi R T)), // NEED / 1000 here Units!!
                //  so
             //      kP := -d(vDot_pair)/dp
             //          = sigma_e*A(T)/delta * y_i * dVspec  >= 0
         //  (x, y, Psat, T and the densities are frozen in the
         //  linearisation; p = p_rgh + rho*g*h, so
                //  d/dp == d/dp_rgh).  In pEqn this goes in as
                        //      + fvm::Sp(vDotP, p_rgh) - vDotP*p_rgh.oldTime()
                        //  making the source backward-Euler in p; summed over
                         //  pairs the realized source is exactly
                        //      SUM_i sigma_e*A_i/delta*dV_i
                        //           *(x_i*Psat_i - y_i*p_new).
                        if (vDotPptr)
                        {
                            volScalarField kP
                            (
                                IOobject
                                (
                                    "kP." + liq.name(),
                                    runTime.timeName(),
                                    mesh_
                                ),
                                accommodationCoeff
                               *yVap
                               *Foam::sqrt
                                (
                                    (Wliq/1000.0)
                                   /(2.0*M_PI*gasconstant*Tsafe)
                                )
                               /int_thickness
                               *max
                                (
                                    1.0/rhoVap - 1.0/rhoLiq,
                                    dimensionedScalar
                                    (
                                        "zeroDV",
                                        dimVolume/dimMass,
                                        0.0
                                    )
                                )
                            );
                            // same gate
                            kP *= gInt;
                            kP.primitiveFieldRef() *= unclampedI;
                            *vDotPptr += kP;
                        }

                        reduce(nClampE, sumOp<label>());
                        reduce(nClampM, sumOp<label>());
                        reduce(nClampV, sumOp<label>());
                        Info<< "    (" << liq.name() << "," << vap.name()
                            << ")  Tboil=" << Tboil.value()
                            << "  max rate magnitude=" << max(mag(r)).value()
                            << " kg/m3/s"
                            << "  clamped cells ENERGY/MASS/VOLUME to maintain boundedness..... probably decrease timestep or pray to the numerical gods: "
                            << nClampE << "/" << nClampM << "/" << nClampV
                            << endl;
                    }
                }
            }
            ++li;
        }
    }


    volScalarField beta // CONDENSED PHASE SUM, for flux limiting step to stay sharp (VOF SHARP)
    (
        IOobject
        (
            "alpha.beta",
            runTime.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        condensate                                   
    );
    beta.oldTime() == beta0_;                        


    const scalar cAlphaBeta
    (
        mesh_.solverDict("alpha").getOrDefault<scalar>("cAlpha", 1.0)
    );

    // un-limited volumetric flux of beta
    surfaceScalarField betaPhiCorr
    (
        "betaPhiCorr",
        fvc::flux(phi_, beta, alphaScheme)
    );

    // compression flux normal to the condensate interface
    {
        surfaceScalarField phic(mag(phi_)/mesh_.magSf());
        surfaceScalarField phir
        (
            min(cAlphaBeta*phic, max(phic))*nHatf(beta, scalar(1) - beta)
        );
        betaPhiCorr += fvc::flux
        (
            -fvc::flux(-phir, (scalar(1) - beta), alpharScheme),
            beta,
            alpharScheme
        );
    }


    volScalarField SpBeta
    (
        IOobject("SpBeta", runTime.timeName(), mesh_),
        mesh_,
        dimensionedScalar("0", dimensionSet(0, 0, -1, 0, 0), 0.0)
    );
    volScalarField SuBeta
    (
        IOobject("SuBeta", runTime.timeName(), mesh_),
        mesh_,
        dimensionedScalar("0", dimensionSet(0, 0, -1, 0, 0), 0.0)
    );
    {
        const scalarField& betaI       = beta.primitiveField();
        const scalarField& dgdtBarI    = dgdtBar.primitiveField();
        const scalarField& dgdtBetaI   = dgdtBetaSum.primitiveField();
        scalarField& SpI = SpBeta.primitiveFieldRef();
        scalarField& SuI = SuBeta.primitiveFieldRef();

        forAll(betaI, celli)
        {
            const scalar b = betaI[celli];

            // dilatation: cancels the beta*div(U) carried in the flux
            SuI[celli] += (-dgdtBarI[celli])*Foam::min(b, scalar(1));

            // relative D = beta*dgdtBar - dgdtBetaSum  (sum of the per-phase
            // relative terms over the condensed group of phases]
            const scalar D = b*dgdtBarI[celli] - dgdtBetaI[celli];
            if (D > 0)
            {
                const scalar w = Foam::max(scalar(1) - b, scalar(1e-4));
                SpI[celli] -= D/w;
                SuI[celli] += D/w;
            }
            else if (D < 0)
            {
                const scalar w = Foam::max(b, scalar(1e-4));
                SpI[celli] += D/w;
            }
        }
    }

    {
        label li = 0;
        for (phaseModel& liq : phases_)
        {
            if (!isGas[li])
            {
                const word vapName = liq.name() + "vapour";
                if (name2idx.found(vapName))
                {
                    // SuMass[li] = -r_i  already; convert mass->volume via rho_liq
                    const volScalarField rhoLiq
                    (
                        max
                        (
                            liq.thermo().rho(),
                            dimensionedScalar("rhoFloorPC", dimDensity, 1e-3)
                        )
                    );
                    SuBeta += SuMass[li]/rhoLiq;     // = -r_i/rho_liq_i
                }
            }
            ++li;
        }
    }

    MULES::limit
    (
        rDeltaT,
        geometricOneField(),
        beta,
        phi_,
        betaPhiCorr,
        SpBeta,
        SuBeta,
        oneField(),
        zeroField(),
        true
    );


    surfaceScalarField betaPhi
    (
        "betaPhi",
        betaPhiCorr + upwind<scalar>(mesh_, phi_).flux(beta)
    );

    {
        volScalarField::Internal SpBetaI
            (IOobject("SpBetaI", runTime.timeName(), mesh_), SpBeta);
        volScalarField::Internal SuBetaI
            (IOobject("SuBetaI", runTime.timeName(), mesh_), SuBeta);

        MULES::explicitSolve
        (
            geometricOneField(),
            beta,
            betaPhi,
            SpBetaI,
            SuBetaI
        );
    }
    beta = max(min(beta, scalar(1)), scalar(0));   // safety pin


    surfaceScalarField gammaPhi("gammaPhi", phi_ - betaPhi);

    const dimensionedScalar betaSmall("betaSmall", dimless, 1e-8);

    rhoPhi_ = dimensionedScalar(dimensionSet(1, 0, -1, 0, 0), Zero);

    {
        label i = 0;
        for (phaseModel& ph : phases_)
        {

            volScalarField hostFrac
            (
                IOobject("hostFrac." + ph.name(), runTime.timeName(), mesh_),
                mesh_,
                dimensionedScalar("hostFrac", dimless, 0.0)
            );
            // START-OF-STEP group fraction, NOT the post-update 
            if (isGas[i])
            {
                hostFrac = scalar(1) - condensate;
            }
            else
            {
                hostFrac = condensate;
            }

            const volScalarField psi
            (
                "psi." + ph.name(),
                min
                (
                    max
                    (
                        c_[i]/max(hostFrac, betaSmall),
                        dimensionedScalar("zeroRho", dimDensity, 0.0)
                    ),
                    ph.thermo().rho()
                )
            );

            // species MASS flux from the sharp group flux  [kg/s]
            const surfaceScalarField& groupPhi = isGas[i] ? gammaPhi : betaPhi;
            surfaceScalarField Fi
            (
                "F" + ph.name(),
                upwind<scalar>(mesh_, groupPhi).flux(psi)
            );

           
            //  DIffusion


            {
                const volScalarField& ai = ph;
                const surfaceScalarField rhoif
                (
                    fvc::interpolate(ph.thermo().rho())
                );

                for (phaseModel& ph2 : phases_)
                {
                    if (&ph2 == &ph) continue;                 // skip self

                    const label j = name2idx[ph2.name()];
                    if (isGas[j] != isGas[i]) continue;        // same group only

                    scalarCoeffSymmDTable::const_iterator dAlpha
                    (
                        dAlphas_.find(interfacePair(ph, ph2))
                    );
                    if (dAlpha == dAlphas_.end()) continue;    // listed pairs only

                    const volScalarField& aj = ph2;
                    const dimensionedScalar Dij("Dij", dimdiff_, dAlpha());


                    Fi +=
                        rhoif
                       *(
                          - Dij*mesh_.magSf()
                           *(
                                fvc::interpolate(aj)*fvc::snGrad(ai)
                              - fvc::interpolate(ai)*fvc::snGrad(aj)
                            )
                        );
                }
            }

            // conservative explicit update of c_i 
            volScalarField::Internal Su
            (
                IOobject("Su_c" + ph.name(), runTime.timeName(), mesh_),
                SuMass[i]
            );

            MULES::explicitSolve
            (
                geometricOneField(),
                c_[i],
                Fi,
                zeroField(),
                Su
            );

            rhoPhi_ += Fi;                          // accumulate mixture mass flux
            ++i;
        }
    }


    //        alpha_i = max(c_i, 0) / rho_i
    //
    //    c_i can numerically carry a small like ~ -(numerical tolerance)

 
    const dimensionedScalar rhoMin("rhoMin", dimDensity, 1e-6);
    const dimensionedScalar cZero("cZero", dimensionSet(1, -3, 0, 0, 0), 0.0);

    volScalarField sumAll
    (
        IOobject("sumAlphaRecovered", runTime.timeName(), mesh_),
        mesh_, dimensionedScalar("0", dimless, 0.0)
    );

    {
        label i = 0;
        for (phaseModel& ph : phases_)
        {
            const volScalarField rhoPh
            (
                max(ph.thermo().rho(), rhoMin)       // positive density stopper/minima
            );

   
            ph == max(c_[i], cZero)/rhoPh;
            sumAll += ph;
            ++i;
        }
    }

    // Volume-closure residual BEFORE rescale
    Info<< "sum(alpha) BEFORE rescale: min = " << min(sumAll).value()
        << ", max = " << max(sumAll).value()
        << ", max|dev| = " << max(mag(sumAll - 1.0)).value() << endl;

    const dimensionedScalar sumFloor("sumFloor", dimless, 1e-30);
    for (phaseModel& ph : phases_)
    {
        ph *= 1.0/max(sumAll, sumFloor);             // global Sum(alpha) = 1
        ph.correctBoundaryConditions();
        ph == max(ph, dimensionedScalar("zero", dimless, 0.0));
    }


    //  DIAGNOSTIC STUFFF
    


    // {
    //     Info<< "Ledger negativity (min c_i, integrated debt):" << endl;
    //     label i = 0;
    //     for (const phaseModel& ph : phases_)
    //     {
    //         const scalarField& cI = c_[i].primitiveField();
    //         const scalar cMin = gMin(cI);
    //         const scalar debt =
    //             gSum(Foam::min(cI, scalar(0))*mesh_.V().field());
    //         Info<< "    " << ph.name() << ": min = " << cMin
    //             << " kg/m3, debt = " << debt << " kg" << endl;
    //         ++i;
    //     }
    // }

    Info<< "Conserved partial-mass(c_i = alpha_i*rho_i):" << endl;
    {
        HashTable<bool> done;
        for (const phaseModel& ph : phases_)
        {
            if (done.found(ph.name())) continue;
            const label i = name2idx[ph.name()];
            const word vapName = ph.name() + "vapour";
            const scalar mLiq =
                gSum(c_[i].primitiveField()*mesh_.V().field());
            if (name2idx.found(vapName))
            {
                const label vi = name2idx[vapName];
                const scalar mVap =
                    gSum(c_[vi].primitiveField()*mesh_.V().field());
                Info<< "    " << ph.name() << "+" << vapName
                    << " = " << mLiq + mVap << " kg (liq " << mLiq
                    << ", vap " << mVap << ")" << endl;
                done.insert(vapName, true);
                done.insert(ph.name(), true);
            }
            else
            {
                Info<< "    " << ph.name() << " = " << mLiq << " kg" << endl;
                done.insert(ph.name(), true);
            }
        }
    }


    volScalarField sumAlpha
    (
        IOobject("sumAlpha", runTime.timeName(), mesh_),
        mesh_, dimensionedScalar("0", dimless, 0.0)
    );
    for (const phaseModel& ph : phases_) sumAlpha += ph;
    Info<< "sum(alpha) after rescale: max|dev| = "
        << max(mag(sumAlpha - 1.0)).value() << endl;

    {
        label i = 0;
        for (const phaseModel& ph : phases_)
        {
            Info<< ph.name() << " alpha min/max/avg = "
                << min(ph).value() << ' ' << max(ph).value() << ' '
                << ph.weightedAverage(mesh_.V()).value() << endl;
            ++i;
        }
    }

    calcAlphas();
    return tPCR;
}


// ************************************************************************* //