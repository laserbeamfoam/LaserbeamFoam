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


//  multiphaseMixtureThermo::solveAlphas   --   the heirarchical volume of fluid aproach
//  one sharp condensed field (MULES + compression keeps the interface crisp) and the
//  actual conserved things are the per species PARTIAL MASSES
//
//        c_i = alpha_i * rho_i        units are  [kg/m^3]
//
//  all mass changes happen to c_i thruogh flux form updates and the exactly paired
//  phase change sources, nothing else is alowed to touch them or conservattion
//  breaks. the alphas are just recoverd from the mass-balance at the end
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

    // wipe last steps latent heat sink and presure coupling, rebuilt fresh
    *massdotterm *= 0.0;

    if (vDotPptr)
    {
        *vDotPptr *= 0.0;
    }


    const dimensionedScalar Tmin("Tmin", dimTemperature, 300.0);// T floor for the rate maths only, stops sqrt/exp FPEs
    const volScalarField Tsafe(max(T_, Tmin)); // safe tempreature for all teh kinetics below, doesnt touch the energy eqn

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

    // per step cap on how much volume phase change can make/remove 
    const scalar etaVol
    (
        phasedictionary.getOrDefault<scalar>("phaseChangeVolLimit", 0.1)
    );

    // old symmetric Hertz-Knudsen accomodation coefficient. kept as the falback

    const scalar accommodationCoeff
    (
        phasedictionary.getOrDefault<scalar>("accommodationCoeff", 1.0)
    );

    // split accomodation: evapouration keeps near 1 but condensattion is
    // throttled (0.01-0.1 ish). real condensation onto a cool surface has to
    // difuse through the air sat next to it which slows it right down and we
    // dont resolve any of that, with symmetric sigma=1 the whole melt pool
    // surface hoovers up the vapour and the keyhole fills with air.....
    const scalar sigmaEvap
    (
        phasedictionary.getOrDefault<scalar>
        (
            "accommodationCoeffEvap",
            accommodationCoeff
        )
    );
    const scalar sigmaCond
    (
        phasedictionary.getOrDefault<scalar>
        (
            "accommodationCoeffCond",
            accommodationCoeff
        )
    );

    Info<< "Phase-change accommodation: sigmaEvap = " << sigmaEvap
        << ", sigmaCond = " << sigmaCond << endl;

    const label nPhases = phases_.size();


    //  phase bookeeping: name to index map and a gas flag list so we arent
    //  asking isGaseous() a millon times

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
    // numerically pin to [0,1] (MULES re-bounds)
    condensate = max(min(condensate, scalar(1)), scalar(0));

    
    //  melt state gate
    //  the interface check further down cant tell melt from solid (epsiln1), so
    //  vapour sat next to a cold solid wall condenses at the full 
    //  rate
    //  so scale r by the melt fraction of the metal thats
    //  actualy present
    scalarField gCond(mesh_.nCells(), 1.0);
    {
        const scalar solidCondFactor =
            phasedictionary.lookupOrDefault<scalar>
            (
                "solidCondensationFactor",
                0.0
            );

        volScalarField TsolC
        (
            IOobject("TsolC", runTime.timeName(), mesh_),
            mesh_,
            dimensionedScalar("Ts0", dimTemperature, 0.0)
        );
        volScalarField TliqC
        (
            IOobject("TliqC", runTime.timeName(), mesh_),
            mesh_,
            dimensionedScalar("Tl0", dimTemperature, 0.0)
        );

        for (const phaseModel& ph : phases_)
        {
            if (ph.isGaseous()) continue;

            // same per-phase dicts update.H reads each step; scope-local
            IOdictionary phdict
            (
                IOobject
                (
                    "thermophysicalProperties." + ph.name(),
                    mesh_.time().constant(),
                    mesh_,
                    IOobject::MUST_READ_IF_MODIFIED
                )
            );

            const dimensionedScalar Ts
            (
                "Ts",
                dimTemperature,
                readScalar(phdict.lookup("TSolidus"))
            );
            const dimensionedScalar Tl
            (
                "Tl",
                dimTemperature,
                readScalar(phdict.lookup("TLiquidus"))
            );

            TsolC += ph*Ts;
            TliqC += ph*Tl;
        }

        const volScalarField betaN
        (
            max(condensate, dimensionedScalar("bSmall", dimless, 1e-6))
        );
        TsolC /= betaN;
        TliqC /= betaN;

        const scalarField& TI  = T_.primitiveField();
        const scalarField& TsI = TsolC.primitiveField();
        const scalarField& TlI = TliqC.primitiveField();

        forAll(gCond, celli)
        {
            const scalar band =
                Foam::max(TlI[celli] - TsI[celli], scalar(1));
            const scalar g = Foam::min
            (
                Foam::max((TI[celli] - TsI[celli])/band, scalar(0)),
                scalar(1)
            );
            gCond[celli] = Foam::max(g, solidCondFactor);
        }
    }


    // r is an interface flux smeared over int_thickness so it should only be
    // where liquid and gas actualy coexist. ungated it misbehaves at both
    // extremes
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


    // once per timestep snapshot. every outer corrector re-solves from the
    // same start of step state (semi PIMPLE) so stash c and beta 
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


    //  compressibility bookeeping for the beta transport
    //    dgdtBar      = sum all phases of alpha*dgdt   (mixture dilattion)
    //    dgdtBetaSum  = same but condensed phases only

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

    //  RAOULT-DALTON multicomponent kinetics:
            //   x_i = mole fraction of species i in the liquid (condensed group only)
    //   y_i = mole fraction in the gas including air.. Dalton, the gas side
        //         only feels its own partial presure y_i*p not the total


   
    const volScalarField rCvMix(this->rCv());   // mixture reciprocal Cv

    const volScalarField rhoMixField(this->rho());


    //
    //  molar densities n_i = max(c_i,0)/W_i then mole fractions within each
    //  group. built straight from teh mass-balance

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

                        // Raoult-Dalton composition of this pair (mass-balance-based)
                        const volScalarField& xLiq = moleFrac[li];
                        const volScalarField& yVap = moleFrac[vi];

                        // Clausius-Clapeyron saturation presure
                        const volScalarField Psat
                        (
                            IOobject("Psat." + liq.name(), runTime.timeName(), mesh_),
                            P0*Foam::exp
                            (
                                (((Wliq/1000.0)*Lvap)/(Tboil*gasconstant))
                               *(1.0 - (Tboil/Tsafe))
                            )
                        );

                        // the one signed rate for this pair [kg/m3/s], positive
                        // evapourates negative condenses. driving force is
                        // x*Psat - y*p, liquid side equilibrium partial presure
                        // vs whats actully sat in the gas
                        volScalarField r
                        (
                            IOobject("r" + liq.name(), runTime.timeName(), mesh_),
                            Foam::sqrt
                            (
                                (Wliq/1000.0)/(2.0*M_PI*gasconstant*Tsafe)
                            )
                           *(xLiq*Psat - yVap*p_)/int_thickness
                        );

                        // apply the split sigma by sign. 
                        {
                            scalarField& rSig = r.primitiveFieldRef();
                            forAll(rSig, celli)
                            {
                                rSig[celli] *=
                                    (rSig[celli] > 0 ? sigmaEvap : sigmaCond);
                            }
                        }

                        
                        // interface gate, see up top. still before the clamps
                        r *= gInt;

                        // solid wall condensation gate (built up top), evap
                        // left alone
                        {
                            scalarField& rI0 = r.primitiveFieldRef();
                            label nGate = 0;
                            forAll(rI0, celli)
                            {
                                if (rI0[celli] < 0 && gCond[celli] < 1.0)
                                {
                                    rI0[celli] *= gCond[celli];
                                    ++nGate;
                                }
                            }
                            reduce(nGate, sumOp<label>());
                            if (nGate)
                            {
                                Info<< "    solid-wall condensation gate ("
                                    << liq.name() << "): " << nGate
                                    << " cell(s)" << endl;
                            }
                        }

                        
                        //  composition consistent saturation tempreature, ie
                                //  the T where x*Psat(T) = y*p and this pairs driving
                            //  force is exactly zero. the energy clamp below then
                                    //  brakes towards the SAME fixed point as the kinetics
                            //  which is the whole point. just invert Psat:
                            //      Tsat = Tboil/(1 - ln(y p/(x P0))/K)
                        //  PAPER REF NEEDED
                        
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

                                //  energy limit: dont let latent heat drive
                                //  T past Tsat this step
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

                                //  mass limit: half of whats there per step,
                                //  max. NB the positivity limiter budget below
                                //  RELIES on this 0.5, dont raise it
                                const scalar rmaxEvap =
                                    0.5*Foam::max(cLiqI[celli], 0.0)*rDeltaT;
                                const scalar rmaxCond =
                                    0.5*Foam::max(cVapI[celli], 0.0)*rDeltaT;
                                if (rc >  rmaxEvap) { rc =  rmaxEvap; ++nClampM; } // evap
                                if (rc < -rmaxCond) { rc = -rmaxCond; ++nClampM; } // cond

                                //  volume limit: etaVol of the cell per step.
                                //  this one scales with 1/dt so it must stay
                                //  under maxPhaseChangeFraction (see note at
                                //  the etaVol read)
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

                        // exactly paired sources, the same one r feeds both
                        // sides so pair mass is conserved to machine precission
                        SuMass[li] -= r;     // liquid loses
                        SuMass[vi] += r;     // vapour gains (paired)

                        // volume made/removed goes to the presure eqn via vDot
                        PCR += r*(1.0/rhoVap - 1.0/rhoLiq);

                        // latent-heat of vap/cond goes to TEqn through -mass_dot ..... maybe i should change the sign, but this is right !!!
                        *massdotterm += rCvMix*Lvap*r;

                        
                        //  implicit presure coupling. the explicit r above is
                        //  frozen at start of step  which lags, 

                        //  slope:
                        //     kP = -d(vDot)/dp = sigma*A/delta*y*dVspec >= 0
                            // and d/dp == d/dp_rgh).
                                //  pEqn adds +fvm::Sp(vDotP,p_rgh) - vDotP*p_rgh.oldTime()
                         //  making the source backward Euler in p, the diagonal
                        //  is allways positive so unconditionaly stabilising same trick as interPhaseChangeFoam
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
                                yVap
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
                            
                            kP *= gInt;

                            
                            {
                                scalarField& kPI = kP.primitiveFieldRef();
                                const scalarField& rIg = r.primitiveField();
                                forAll(kPI, celli)
                                {
                                    if (rIg[celli] < 0)
                                    {
                                        kPI[celli] *= sigmaCond*gCond[celli];
                                    }
                                    else
                                    {
                                        kPI[celli] *= sigmaEvap;
                                    }
                                }
                            }
                            kP.primitiveFieldRef() *= unclampedI;
                            *vDotPptr += kP;
                        }

                        reduce(nClampE, sumOp<label>());
                        reduce(nClampM, sumOp<label>());
                        reduce(nClampV, sumOp<label>());

                        
                        const scalarField& rBud = r.primitiveField();
                        const scalar evapKgPerS =
                            gSum(pos(rBud)*rBud*mesh_.V().field());
                        const scalar condKgPerS =
                            gSum(neg(rBud)*rBud*mesh_.V().field());

                        Info<< "    (" << liq.name() << "," << vap.name()
                            << ")  Tboil=" << Tboil.value()
                            << "  max rate magnitude=" << max(mag(r)).value()
                            << " kg/m3/s"
                            << "  clamped cells ENERGY/MASS/VOLUME to maintain boundedness..... probably decrease timestep or pray to the numerical gods: "
                            << nClampE << "/" << nClampM << "/" << nClampV
                            << endl;
                        Info<< "        integrated evap = " << evapKgPerS
                            << " kg/s, cond = " << condKgPerS
                            << " kg/s (cond is <= 0)" << endl;
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

    // raw un-limited beta flux
    surfaceScalarField betaPhiCorr
    (
        "betaPhiCorr",
        fvc::flux(phi_, beta, alphaScheme)
    );

    // interFoam style compresion flux normal to the interface, keeps it sharp
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

            // dilation bit: cancels the beta*div(U) hiding in the flux
            SuI[celli] += (-dgdtBarI[celli])*Foam::min(b, scalar(1));

            // relative compressibility between the two groups, pushed to
            // whichever side keeps beta bounded... this is nice if i do say so myself
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
                    // SuMass already holds -r for the liquid, just convert
                    // mass to volume with rho_liq 
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
    beta = max(min(beta, scalar(1)), scalar(0));   // safety 


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
            // start of step group fraction
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

            // species mass flux 
            const surfaceScalarField& groupPhi = isGas[i] ? gammaPhi : betaPhi;
            surfaceScalarField Fi
            (
                "F" + ph.name(),
                upwind<scalar>(mesh_, groupPhi).flux(psi)
            );

           
       

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
                    if (isGas[j] != isGas[i]) continue;        

                    scalarCoeffSymmDTable::const_iterator dAlpha
                    (
                        dAlphas_.find(interfacePair(ph, ph2))
                    );
                    if (dAlpha == dAlphas_.end()) continue;    

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

            
            // donor cell positivity limiter on the TOTAL species flux
            // (advection + diffusion) MULEs bounds condensed (beta)
            
            {
                const scalarField& cI  = c_[i].primitiveField();
                const scalarField& SuI = SuMass[i].primitiveField();
                const scalarField& V   = mesh_.V().field();

                scalarField outflow(mesh_.nCells(), 0.0);

                const labelUList& own = mesh_.owner();
                const labelUList& nei = mesh_.neighbour();

                {
                    const scalarField& FiIf = Fi.primitiveField();
                    forAll(FiIf, facei)
                    {
                        if (FiIf[facei] > 0)
                        {
                            outflow[own[facei]] += FiIf[facei];
                        }
                        else
                        {
                            outflow[nei[facei]] -= FiIf[facei];
                        }
                    }
                }

                forAll(Fi.boundaryField(), patchi)
                {
                    const fvsPatchScalarField& pF =
                        Fi.boundaryField()[patchi];
                    const labelUList& fc = pF.patch().faceCells();
                    forAll(pF, bfacei)
                    {
                        if (pF[bfacei] > 0)
                        {
                            outflow[fc[bfacei]] += pF[bfacei];
                        }
                    }
                }


                volScalarField lam
                (
                    IOobject
                    (
                        "lambda." + ph.name(),
                        runTime.timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar("one", dimless, 1.0)
                );
                scalarField& lamI = lam.primitiveFieldRef();

                label nLim = 0;
                scalar lamMin = 1.0;
                forAll(lamI, celli)
                {
                    if (outflow[celli] > VSMALL)
                    {
                        const scalar budget =
                            Foam::max(cI[celli], scalar(0))
                           *V[celli]*rDeltaT
                          + Foam::min(SuI[celli], scalar(0))*V[celli];

                        const scalar l = Foam::min
                        (
                            scalar(1),
                            Foam::max(budget, scalar(0))/outflow[celli]
                        );

                        if (l < 1.0)
                        {
                            ++nLim;
                            lamMin = Foam::min(lamMin, l);
                        }
                        lamI[celli] = l;
                    }
                }

                lam.correctBoundaryConditions();   // fill processor halos

                {
                    scalarField& FiI = Fi.primitiveFieldRef();
                    forAll(FiI, facei)
                    {
                        FiI[facei] *=
                            (FiI[facei] > 0)
                          ? lamI[own[facei]]
                          : lamI[nei[facei]];
                    }
                }

                auto& FiBf = Fi.boundaryFieldRef();
                forAll(FiBf, patchi)
                {
                    fvsPatchScalarField& pF = FiBf[patchi];
                    const labelUList& fc = pF.patch().faceCells();

                    if (pF.patch().coupled())
                    {
                        const scalarField lamNei
                        (
                            lam.boundaryField()[patchi]
                           .patchNeighbourField()
                        );
                        forAll(pF, bfacei)
                        {
                            pF[bfacei] *=
                                (pF[bfacei] > 0)
                              ? lamI[fc[bfacei]]
                              : lamNei[bfacei];
                        }
                    }
                    else
                    {
                        forAll(pF, bfacei)
                        {
                            if (pF[bfacei] > 0)
                            {
                                pF[bfacei] *= lamI[fc[bfacei]];
                            }
                        }
                    }
                }

                reduce(nLim, sumOp<label>());
                reduce(lamMin, minOp<scalar>());
                if (nLim)
                {
                    Info<< "positivity limiter (" << ph.name() << "): "
                        << nLim << " cell(s), min lambda = "
                        << lamMin << endl;
                }
            }

            // THE conservative update of c_i, flux plus the paired sources
            // and nothing else ever
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


    //  now recover alphas from the mass-balance. c_i can carry a tiny negative
    //  from roundoff hence the max(..,0) everywhere

 

        //  the 0.5 condensation clamp is multiplicative so vapour c
            // decays exponentialy and never actually reaches zero, every cell a
    // vapour ever visited keeps a meaningles residue 
        //   snap |c| below cPurgeTol to exactly 0. NB this is the ONE deliberate
    // conservation leak in the whole scheme, its bounded and the purged mass
    // gets printed every step.. cPurgeTol 0 turns it off, but its absolutely minisucle 
    {
        const scalar cPurgeTol =
            phasedictionary.getOrDefault<scalar>("cPurgeTol", 1e-10);

        if (cPurgeTol > 0)
        {
            label nPurge = 0;
            scalar purgedMass = 0.0;
            const scalarField& V = mesh_.V().field();

            forAll(c_, i)
            {
                scalarField& cI = c_[i].primitiveFieldRef();
                forAll(cI, celli)
                {
                    const scalar cv = cI[celli];
                    if (cv != 0 && mag(cv) < cPurgeTol)
                    {
                        purgedMass += mag(cv)*V[celli];
                        cI[celli] = 0.0;
                        ++nPurge;
                    }
                }
            }

            reduce(nPurge, sumOp<label>());
            reduce(purgedMass, sumOp<scalar>());
            if (nPurge)
            {
                Info<< "mass-balance purge (|c| < " << cPurgeTol << " kg/m3): "
                    << nPurge << " entries, " << purgedMass
                    << " kg absolute" << endl;
            }
        }
    }

    const dimensionedScalar cZero("cZero", dimensionSet(1, -3, 0, 0, 0), 0.0);


    //  HIERARCHICAL VOLUME-CLOSURE RECOVERY  (global rescale was shit))
    //
    //  alpha = c/rho only sums to 1 if the transported masses match the EOS
    //  densities at the current p and T, and they never quite do (species
    //  move once per step with last steps fluxes while p moves by bars at the
    //  spot) so a defect builds every step
    // 
        //  condensed alphas come straight from the mass-balance 

    const dimensionedScalar rhoMinRec("rhoMinRec", dimDensity, 1e-3);

    // closure trust thresholds
    //   closureGasFracMin: the gas mass-balance must account for at least this
    //     fraction of the volume  asked to fill 
    //   closureCondNearFull: above this much condensed fill the leftover is
    //     splitting defect not a real void, condensed just absorbs it 
    const scalar gasFracMin =
        phasedictionary.getOrDefault<scalar>("closureGasFracMin", 0.01);
    const scalar condNearFull =
        phasedictionary.getOrDefault<scalar>("closureCondNearFull", 0.95);

    volScalarField vCond
    (
        IOobject("vCondImplied", runTime.timeName(), mesh_),
        mesh_, dimensionedScalar("0", dimless, 0.0)
    );
    volScalarField vGas
    (
        IOobject("vGasImplied", runTime.timeName(), mesh_),
        mesh_, dimensionedScalar("0", dimless, 0.0)
    );

    PtrList<volScalarField> vImp(nPhases);
    {
        label i = 0;
        for (const phaseModel& ph : phases_)
        {
            const volScalarField rhoPh(max(ph.thermo().rho(), rhoMinRec));

            vImp.set
            (
                i,
                new volScalarField
                (
                    "vImp." + ph.name(),
                    max(c_[i], cZero)/rhoPh
                )
            );

            if (isGas[i]) vGas += vImp[i];
            else          vCond += vImp[i];
            ++i;
        }
    }

    const volScalarField sumImplied(vCond + vGas);

    // mass-balance implied volume vs actual cell volume, internal field only
    Info<< "volume closure defect (mass-balance vs EOS): min = "
        << gMin(sumImplied.primitiveField())
        << ", max = " << gMax(sumImplied.primitiveField())
        << ", max|dev| = "
        << gMax(mag(sumImplied.primitiveField() - 1.0)()) << endl;

    {
        // previous gas fractions.... composition-memory for fallback
        label nGasPh = 0;
        forAll(isGas, i)
        {
            if (isGas[i]) ++nGasPh;
        }

        PtrList<scalarField> oldGas(nPhases);
        scalarField oldGasSum(mesh_.nCells(), 0.0);
        {
            label i = 0;
            for (const phaseModel& ph : phases_)
            {
                if (isGas[i])
                {
                    oldGas.set(i, new scalarField(ph.primitiveField()));
                    oldGasSum += oldGas[i];
                }
                ++i;
            }
        }

        const scalarField& vCondI = vCond.primitiveField();
        const scalarField& vGasI  = vGas.primitiveField();

        // condensed group: mass-balance volumes
        {
            label i = 0;
            for (phaseModel& ph : phases_)
            {
                if (!isGas[i])
                {
                    scalarField& aI = ph.primitiveFieldRef();
                    const scalarField& vI = vImp[i].primitiveField();
                    forAll(aI, celli)
                    {
                        const scalar vc = vCondI[celli];
                        const scalar Vg = 1.0 - Foam::min(vc, scalar(1));
                        const scalar vg = vGasI[celli];

                        scalar s = 1.0;
                        if (vc > 1.0)
                        {
                            s = 1.0/vc;
                        }
                        else if
                        (
                            Vg > 0
                         && vg < gasFracMin*Vg
                         && vc >= condNearFull
                        )
                        {
                            // absorb the splitting defect; the gas keeps
                            // its contribution vg, so the cell
                            // still closes exactly: s*vc + vg = 1
                            s = (1.0 - vg)/vc;
                        }

                        aI[celli] = s*vI[celli];
                    }
                }
                ++i;
            }
        }

        // gas group: fills the remaining volume in mass-balance proportion 
        label nVoidFill = 0;
        {
            label i = 0;
            for (phaseModel& ph : phases_)
            {
                if (isGas[i])
                {
                    scalarField& aI = ph.primitiveFieldRef();
                    const scalarField& vI = vImp[i].primitiveField();
                    const scalarField& gI = oldGas[i];

                    forAll(aI, celli)
                    {
                        const scalar vc = vCondI[celli];
                        const scalar Vg = 1.0 - Foam::min(vc, scalar(1));
                        const scalar vg = vGasI[celli];

                        if (Vg <= 0)
                        {
                            aI[celli] = 0.0;
                        }
                        else if (vg >= gasFracMin*Vg)
                        {
                            // trusted mass-balance composition
                            aI[celli] = Vg*vI[celli]/vg;
                        }
                        else if (vc >= condNearFull)
                        {
                            // condensed absorbed the defect above
                            aI[celli] = vI[celli];
                        }
                        else if (oldGasSum[celli] > VSMALL)
                        {
                            aI[celli] = Vg*gI[celli]/oldGasSum[celli];
                            if (Vg > SMALL) ++nVoidFill;
                        }
                        else
                        {
                            aI[celli] = Vg/nGasPh;
                            if (Vg > SMALL) ++nVoidFill;
                        }
                    }
                }
                ++i;
            }
        }

        label nOverfullCond = 0;
        label nAbsorb = 0;
        forAll(vCondI, celli)
        {
            const scalar vc = vCondI[celli];
            if (vc > 1.0)
            {
                ++nOverfullCond;
            }
            else
            {
                const scalar Vg = 1.0 - vc;
                if
                (
                    Vg > 0
                 && vGasI[celli] < gasFracMin*Vg
                 && vc >= condNearFull
                )
                {
                    ++nAbsorb;
                }
            }
        }

        for (phaseModel& ph : phases_)
        {
            ph.correctBoundaryConditions();
        }

        reduce(nOverfullCond, sumOp<label>());
        reduce(nAbsorb, sumOp<label>());
        reduce(nVoidFill, sumOp<label>());
        if (nOverfullCond || nAbsorb || nVoidFill)
        {
            Info<< "hierarchical recovery: condensed-overfull cells = "
                << nOverfullCond << ", defect-absorbed cells = " << nAbsorb
                << ", gas-void fallback cells = "
                << nVoidFill/max(nGasPh, label(1)) << endl;
        }
    }

    // dump the closure diff into the pressuqre equation term. 
    {
        const scalar closureRelax =
            phasedictionary.getOrDefault<scalar>("closureRelax", 0.5);
        const scalar closureVolLimit =
            phasedictionary.getOrDefault<scalar>("closureVolLimit", 0.02);

        if (closureRelax > 0)
        {
            const scalarField& sI = sumImplied.primitiveField();
            const scalar srcMax = closureVolLimit*rDeltaT;

            scalarField src(mesh_.nCells(), 0.0);
            scalarField unclampedC(mesh_.nCells(), 1.0);
            label nClampC = 0;

            forAll(src, celli)
            {
                scalar s = closureRelax*(sI[celli] - 1.0)*rDeltaT;
                if (s >  srcMax) { s =  srcMax; unclampedC[celli] = 0; ++nClampC; }
                if (s < -srcMax) { s = -srcMax; unclampedC[celli] = 0; ++nClampC; }
                src[celli] = s;
            }

            PCR.primitiveFieldRef() += src;

            if (vDotPptr)
            {
                scalarField kPc(mesh_.nCells(), 0.0);
                label i = 0;
                for (const phaseModel& ph : phases_)
                {
                    const volScalarField rhoPh
                    (
                        max(ph.thermo().rho(), rhoMinRec)
                    );
                    const scalarField& rI   = rhoPh.primitiveField();
                    const scalarField& psiI =
                        ph.thermo().psi().primitiveField();
                    const scalarField& cI   = c_[i].primitiveField();

                    forAll(kPc, celli)
                    {
                        kPc[celli] +=
                            Foam::max(cI[celli], scalar(0))*psiI[celli]
                           /Foam::sqr(rI[celli]);
                    }
                    ++i;
                }

                vDotPptr->primitiveFieldRef() +=
                    closureRelax*rDeltaT*kPc*unclampedC;
            }

            reduce(nClampC, sumOp<label>());
            Info<< "closure feedback: relax = " << closureRelax
                << ", capped cells = " << nClampC << endl;
        }
    }


    //  DIAGNOSTIC STUFFF
    


    // {
    //     Info<< "mass-balance negativity (min c_i, integrated debt):" << endl;
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


    {
        scalarField sumChk(mesh_.nCells(), 0.0);
        for (const phaseModel& ph : phases_)
        {
            sumChk += ph.primitiveField();
        }
        // should be machine zero by construction, anything else means the
        // closure maths is broke
        Info<< "sum(alpha) after hierarchical recovery: max|dev| = "
            << gMax(mag(sumChk - 1.0)()) << endl;
    }

    {
        label i = 0;
        for (const phaseModel& ph : phases_)
        {
            Info<< ph.name() << " alpha min/max/avg = "
                << gMin(ph.primitiveField()) << ' '
                << gMax(ph.primitiveField()) << ' '
                << ph.weightedAverage(mesh_.V()).value() << endl;
            ++i;
        }
    }

    // rebuild the stacked alphas indicator for paraview
    calcAlphas();
    return tPCR;
}


// ************************************************************************* //