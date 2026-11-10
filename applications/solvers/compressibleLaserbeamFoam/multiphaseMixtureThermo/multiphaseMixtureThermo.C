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
#include "fvcSurfaceIntegrate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDdt.H"
#include "fvcFlux.H"
#include "fvcMeshPhi.H"
#include "surfaceInterpolate.H"
#include "unitConversion.H"
#include "regionSplit.H"
#include <cmath>

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

    rhoCvPhi_
    (
        IOobject
        (
            "rhoCvPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimEnergy/dimTemperature/dimTime, Zero)
    ),

    rhoCvTPhi_
    (
        IOobject
        (
            "rhoCvTPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimEnergy/dimTime, Zero)
    ),

    divPhiLedger_
    (
        IOobject
        (
            "divPhiLedger",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    ),

    pcMomentum_
    (
        IOobject
        (
            "phaseChangeMomentum",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(dimForce/dimVolume, Zero)
    ),

    pcGasDil_
    (
        IOobject
        (
            "phaseChangeGasDilatation",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    ),

    pcCondDil_
    (
        IOobject
        (
            "phaseChangeCondensedDilatation",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
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
    rhoCvPhi_.setOriented();
    rhoCvTPhi_.setOriented();
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

    
    if (!(mesh_.time().timeIndex() > 0))
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

        rhoPhi_oldTime.write();

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


void Foam::multiphaseMixtureThermo::correctRho(const volScalarField& dp)
{
    for (phaseModel& phase : phases_)
    {
        volScalarField& rhoPh = phase.thermo().rho();
        rhoPh += phase.thermo().psi()*dp;
        // was 1e-4.. 1/rho feeds the p matrix coefficients and the extra
        // decade of conditioning bought nothing, match the other floors
        rhoPh = max(rhoPh, dimensionedScalar("rhoMinEOS", dimDensity, 1e-3));
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


const Foam::volScalarField& Foam::multiphaseMixtureThermo::c
(
    const word& phaseName
) const
{
    label i = 0;
    for (const phaseModel& ph : phases_)
    {
        if (ph.name() == phaseName)
        {
            return c_[i];
        }
        ++i;
    }

    FatalErrorInFunction
        << "No phase named " << phaseName << " in the mixture"
        << exit(FatalError);

    return c_[0];
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::rhoM() const
{
    // not named "rhoM": the solver owns a registered field of that name
    auto trhoM = volScalarField::New
    (
        "rhoM.ledger",
        IOobject::NO_REGISTER,
        mesh_,
        dimensionedScalar(dimDensity, Zero)
    );
    volScalarField& r = trhoM.ref();

    for (const volScalarField& ci : c_)
    {
        r += ci;
    }

    return trhoM;
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

            // CSF, Brackbill et al. 1992, https://doi.org/10.1016/0021-9991(92)90240-Y (sigma constant)
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
    volScalarField* vDotDtPtr,
    volScalarField* closureSrcPtr
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

    PCR = solveAlphas(massdotterm, vDotDtPtr, closureSrcPtr);

    writeOldTimeValues();

    return tPCR;
}


Foam::tmp<Foam::surfaceVectorField> Foam::multiphaseMixtureThermo::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
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


// solveAlphas: one sharp condensed field beta (MULES + compression) carries the interface; the conserved
// quantities are the partial masses c_i = alpha_i*rho_i [kg/m^3], changed only by flux-form updates and
// exactly paired phase-change sources

Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::solveAlphas
(
    volScalarField* massdotterm,
    volScalarField* vDotDtPtr,
    volScalarField* closureSrcPtr
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

    // wipe last steps latent heat sink, momentum transfer and dt-control dilatation, rebuilt fresh
    *massdotterm *= 0.0;
    pcMomentum_ = dimensionedVector(pcMomentum_.dimensions(), Zero);
    pcGasDil_ = dimensionedScalar(pcGasDil_.dimensions(), Zero);
    pcCondDil_ = dimensionedScalar(pcCondDil_.dimensions(), Zero);

    if (vDotDtPtr)
    {
        *vDotDtPtr *= 0.0;
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
    // (interface_thickness, localInterfaceThickness and interfaceCellFactor are no
    // longer read: the surface flux is placed by the interface area density below)

    // per-step cap on the volume phase change makes or removes: twice the dt control's target (setDeltaT.H)
    const scalar etaVol
    (
        phasedictionary.getOrDefault<scalar>
        (
            "phaseChangeVolLimit",
            2*runTime.controlDict().getOrDefault<scalar>("maxPhaseChangeFraction", 0.25)
        )
    );

    // old symmetric Hertz-Knudsen accomodation coefficient, kept as the falback
    // default so ancient dicts still run. dont use it symmetric though, see below

    const scalar accommodationCoeff
    (
        phasedictionary.getOrDefault<scalar>("accommodationCoeff", 1.0)
    );

    // split accomodation
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

    
    //  melt state gate for the condensing branch  the interface check further down cant tell melt from solid, so vapour sat next to a cold solid wall condenses at the full kinetic rate and the wall  turns into a vaccum pump that inhales air down the crater
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


    // Phase change: the Hertz-Knudsen flux j acts on the anchors' condensed phase through the area density
    // aI = 2*beta*|grad(beta)| (unit integral across the surface). The liquid goes from the top; the vapour
    // is laid where there is gas to push (kPC) and drawn from there on condensation (tests 07-09)
    const volScalarField magGradBeta(mag(fvc::grad(condensate)));
    volScalarField aI
    (
        IOobject("phaseChangeArea", runTime.timeName(), mesh_),
        2.0*condensate*magGradBeta
    );

    // f seen from the interface: the w-weighted mean over the cell's faces of the two-cell face means;
    // a cell without the phase takes its neighbours' value
    auto interfaceState = [&](const volScalarField& f, const volScalarField& w)
    {
        const dimensionedScalar smallA("smallA", dimArea, VSMALL);
        return tmp<volScalarField>::New
        (
            "interfaceState(" + f.name() + ")",
            fvc::surfaceSum(fvc::interpolate(w*f)*mesh_.magSf())
           /(fvc::surfaceSum(fvc::interpolate(w)*mesh_.magSf()) + smallA)
        );
    };

    // no phase change in cells touching phaseChangeMaskPatches (open truncations of the atmosphere): vapour
    // that would condense far downstream condensed AT the patch and froze a crust across the exhaust, so it
    // leaves as vapour. Masked cells neither anchor nor receive
    volScalarField pcOn
    (
        IOobject("phaseChangeOn", runTime.timeName(), mesh_),
        mesh_,
        dimensionedScalar("pcOn", dimless, 1.0)
    );
    {
        wordList pcMaskNames
        (
            phasedictionary.lookupOrDefault<wordList>
            (
                "phaseChangeMaskPatches",
                wordList()
            )
        );

        label nPCMask = 0;
        scalarField& onI = pcOn.primitiveFieldRef();
        const scalarField& aII = aI.primitiveField();

        for (const word& pn : pcMaskNames)
        {
            const label patchi = mesh_.boundaryMesh().findPatchID(pn);
            if (patchi < 0)
            {
                WarningInFunction
                    << "phaseChangeMaskPatches: patch " << pn
                    << " not found, skipping" << endl;
                continue;
            }
            const labelUList& fc = mesh_.boundary()[patchi].faceCells();
            forAll(fc, bfacei)
            {
                const label celli = fc[bfacei];
                if (onI[celli] > 0 && aII[celli] > 0) ++nPCMask;
                onI[celli] = 0.0;
            }
        }
        reduce(nPCMask, sumOp<label>());
        if (nPCMask)
        {
            Info<< "phase change masked in " << nPCMask
                << " boundary cell(s)" << endl;
        }
    }
    pcOn.correctBoundaryConditions();

    // face shares [m2]: the drop in beta across the face times its area, from more to less condensed (owner to
    // neighbour positive, outward on coupled patches); WPC sums a cell's outgoing shares. The drop comes from
    // the face mean of the Gauss gradients, which is blind to odd-even noise in beta (test 15)
    const volVectorField gradBetaL
    (
        "gradBetaL",
        fvc::surfaceIntegrate(mesh_.Sf()*linearInterpolate(condensate))
    );
    const surfaceScalarField dropL
    (
        "dropL",
        -(linearInterpolate(gradBetaL) & mesh_.delta())
    );
    surfaceScalarField wPC
    (
        IOobject("phaseChangeShare", runTime.timeName(), mesh_),
        mesh_,
        dimensionedScalar("wPC", dimArea, Zero)
    );
    volScalarField WPC
    (
        IOobject("phaseChangeShareSum", runTime.timeName(), mesh_),
        mesh_,
        dimensionedScalar("WPC", dimArea, Zero)
    );
    {
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        const scalarField& bI = condensate.primitiveField();
        const scalarField& onI = pcOn.primitiveField();
        const scalarField& magSfI = mesh_.magSf().primitiveField();
        scalarField& wI = wPC.primitiveFieldRef();
        scalarField& WI = WPC.primitiveFieldRef();

        const scalarField& dLI = dropL.primitiveField();
        forAll(own, facei)
        {
            const label o = own[facei];
            const label n = nei[facei];
            const scalar d = bI[o] - bI[n];
            if (d > 0)
            {
                wI[facei] = Foam::max(dLI[facei], 0.0)*magSfI[facei]*onI[n];
                WI[o] += wI[facei];
            }
            else if (d < 0)
            {
                wI[facei] = -Foam::max(-dLI[facei], 0.0)*magSfI[facei]*onI[o];
                WI[n] -= wI[facei];
            }
        }

        forAll(wPC.boundaryField(), patchi)
        {
            const fvPatch& fvp = mesh_.boundary()[patchi];
            scalarField& wP = wPC.boundaryFieldRef()[patchi];
            wP = 0.0;
            if (!fvp.coupled()) continue;

            const labelUList& fc = fvp.faceCells();
            const scalarField bNb
            (
                condensate.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField onNb
            (
                pcOn.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField& mSf = fvp.magSf();
            const scalarField& dLP = dropL.boundaryField()[patchi];
            forAll(fc, i)
            {
                const scalar d = bI[fc[i]] - bNb[i];
                if (d > 0)
                {
                    wP[i] = Foam::max(dLP[i], 0.0)*mSf[i]*onNb[i];
                    WI[fc[i]] += wP[i];
                }
                else if (d < 0)
                {
                    wP[i] = -Foam::max(-dLP[i], 0.0)*mSf[i]*onI[fc[i]];
                }
            }
        }

        // an anchor with nowhere to lay the vapour has no phase change
        scalarField& aII = aI.primitiveFieldRef();
        forAll(aII, celli)
        {
            if (!(WI[celli] > 0) || !(onI[celli] > 0)) aII[celli] = 0.0;
        }
    }

    // what each cell receives through the faces where it has less condensed phase:
    // the sender's value per unit share times the share, summed
    auto pcReceive = [&](const volScalarField& q)
    {
        tmp<scalarField> tR(new scalarField(mesh_.nCells(), 0.0));
        scalarField& R = tR.ref();
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        const scalarField& wI = wPC.primitiveField();
        const scalarField& qI = q.primitiveField();
        forAll(own, facei)
        {
            const scalar w = wI[facei];
            if (w > 0) R[nei[facei]] += qI[own[facei]]*w;
            else if (w < 0) R[own[facei]] -= qI[nei[facei]]*w;
        }
        forAll(wPC.boundaryField(), patchi)
        {
            if (!mesh_.boundary()[patchi].coupled()) continue;
            const scalarField& wP = wPC.boundaryField()[patchi];
            const labelUList& fc = mesh_.boundary()[patchi].faceCells();
            const scalarField qNb(q.boundaryField()[patchi].patchNeighbourField());
            forAll(fc, i)
            {
                if (wP[i] < 0) R[fc[i]] -= qNb[i]*wP[i];
            }
        }
        return tR;
    };

    // the mirror: for each sender, the share-weighted sum of a receiver value
    auto pcSendWeighted = [&](const volScalarField& s)
    {
        tmp<scalarField> tS(new scalarField(mesh_.nCells(), 0.0));
        scalarField& S = tS.ref();
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        const scalarField& wI = wPC.primitiveField();
        const scalarField& sI = s.primitiveField();
        forAll(own, facei)
        {
            const scalar w = wI[facei];
            if (w > 0) S[own[facei]] += w*sI[nei[facei]];
            else if (w < 0) S[nei[facei]] -= w*sI[own[facei]];
        }
        forAll(wPC.boundaryField(), patchi)
        {
            if (!mesh_.boundary()[patchi].coupled()) continue;
            const scalarField& wP = wPC.boundaryField()[patchi];
            const labelUList& fc = mesh_.boundary()[patchi].faceCells();
            const scalarField sNb(s.boundaryField()[patchi].patchNeighbourField());
            forAll(fc, i)
            {
                if (wP[i] > 0) S[fc[i]] += wP[i]*sNb[i];
            }
        }
        return tS;
    };


    // kPC: the share of its vapour a cell keeps (the rest passes on through its shares): the gas share of the
    // momentum under Fauske's critical-flow slip u_g/u_l = sqrt(rho_l/rho_g), ANL-6633 (1962),
    // https://digital.library.unt.edu/ark:/67531/metadc1035802/m2/1/high_res_d/4749073.pdf. Volume or mass shares alone failed tests 14-16
    volScalarField kPC
    (
        IOobject("phaseChangeKeep", runTime.timeName(), mesh_),
        mesh_,
        dimensionedScalar("kPC", dimless, Zero)
    );
    {
        scalarField sG(mesh_.nCells(), 0.0);
        scalarField sA(mesh_.nCells(), 0.0);
        label i = 0;
        for (const phaseModel& ph : phases_)
        {
            // alpha*sqrt(rho) = sqrt(alpha*c) from the ledger
            const scalarField s
            (
                Foam::sqrt
                (
                    Foam::max(ph.primitiveField(), 0.0)
                   *Foam::max(c_[i].primitiveField(), 0.0)
                )
            );
            sA += s;
            if (isGas[i]) sG += s;
            ++i;
        }
        scalarField& kI = kPC.primitiveFieldRef();
        forAll(kI, celli)
        {
            kI[celli] = (sA[celli] > 0) ? sG[celli]/sA[celli] : 1.0;
        }
        kPC.correctBoundaryConditions();
    }

    // the band in decreasing beta: shares only run downhill, so every sender comes before
    // its receivers
    labelList pcOrder;
    {
        const scalarField& WI = WPC.primitiveField();
        const scalarField& bI = condensate.primitiveField();
        DynamicList<label> band;
        forAll(WI, celli)
        {
            if (WI[celli] > 0) band.append(celli);
        }
        scalarField key(band.size());
        forAll(band, i) key[i] = -bI[band[i]];
        labelList idx;
        sortedOrder(key, idx);
        pcOrder.setSize(band.size());
        forAll(idx, i) pcOrder[i] = band[idx[i]];
    }
    bool pcCoupled = false;
    forAll(mesh_.boundary(), patchi)
    {
        if (mesh_.boundary()[patchi].coupled()) pcCoupled = true;
    }
    reduce(pcCoupled, orOp<bool>());

    // down the shares: Y = own + what the senders pass on, (1 - k)*Y*w/W through each
    // share. Across coupled faces the inflow is iterated until it settles
    auto pcSweepDown = [&](const List<scalarField>& own)
    {
        const label nq = own.size();
        const labelUList& ow = mesh_.owner();
        const labelUList& ne = mesh_.neighbour();
        const cellList& cells = mesh_.cells();
        const label nInt = mesh_.nInternalFaces();
        const scalarField& wI = wPC.primitiveField();
        const scalarField& WI = WPC.primitiveField();
        const scalarField& kI = kPC.primitiveField();
        List<scalarField> Y(own);
        PtrList<volScalarField> g(nq);
        for (label q = 0; q < nq; ++q)
        {
            g.set
            (
                q,
                new volScalarField
                (
                    IOobject("pcG" + Foam::name(q), runTime.timeName(), mesh_),
                    mesh_,
                    dimensionedScalar(dimless, Zero)
                )
            );
        }
        scalar tol = 0;
        for (label q = 0; q < nq; ++q) tol = Foam::max(tol, gSum(mag(own[q])));
        tol = 1e-12*tol + VSMALL;
        label sweep = 0;
        scalar change = GREAT;
        for (; sweep < 64 && change > tol; ++sweep)
        {
            List<scalarField> I(nq, scalarField(mesh_.nCells(), 0.0));
            List<List<scalarField>> gUsed(nq);
            for (label q = 0; q < nq; ++q)
            {
                gUsed[q].setSize(mesh_.boundary().size());
                forAll(wPC.boundaryField(), patchi)
                {
                    if (!mesh_.boundary()[patchi].coupled()) continue;
                    const scalarField& wP = wPC.boundaryField()[patchi];
                    const labelUList& fc = mesh_.boundary()[patchi].faceCells();
                    gUsed[q][patchi] = g[q].boundaryField()[patchi].patchNeighbourField();
                    forAll(fc, i)
                    {
                        if (wP[i] < 0) I[q][fc[i]] -= gUsed[q][patchi][i]*wP[i];
                    }
                }
            }
            for (const label c : pcOrder)
            {
                const scalar pass = (1.0 - kI[c])/WI[c];
                for (label q = 0; q < nq; ++q)
                {
                    Y[q][c] = own[q][c] + I[q][c];
                    g[q].primitiveFieldRef()[c] = pass*Y[q][c];
                }
                for (const label facei : cells[c])
                {
                    if (facei >= nInt) continue;
                    const scalar w = wI[facei];
                    label n = -1;
                    if (w > 0 && ow[facei] == c) n = ne[facei];
                    else if (w < 0 && ne[facei] == c) n = ow[facei];
                    if (n < 0) continue;
                    for (label q = 0; q < nq; ++q)
                    {
                        I[q][n] += g[q].primitiveField()[c]*mag(w);
                    }
                }
            }
            forAll(WI, celli)
            {
                if (WI[celli] > 0) continue;
                for (label q = 0; q < nq; ++q) Y[q][celli] = own[q][celli] + I[q][celli];
            }
            if (!pcCoupled) break;
            change = 0;
            for (label q = 0; q < nq; ++q)
            {
                g[q].correctBoundaryConditions();
                forAll(wPC.boundaryField(), patchi)
                {
                    if (!mesh_.boundary()[patchi].coupled()) continue;
                    const scalarField& wP = wPC.boundaryField()[patchi];
                    const scalarField gNb(g[q].boundaryField()[patchi].patchNeighbourField());
                    forAll(wP, i)
                    {
                        change = Foam::max(change, mag((gNb[i] - gUsed[q][patchi][i])*wP[i]));
                    }
                }
            }
            reduce(change, maxOp<scalar>());
        }
        if (pcCoupled && change > tol)
        {
            Info<< "phase change: coupled sweep (down) not converged after " << sweep
                << " sweeps, change " << change << endl;
        }
        return Y;
    };

    // up the shares, the mirror: X = b + (1 - k)*sum over the cell's shares of (w/W)*X of
    // the receiver. With b = k it is the part of a cell's vapour that is laid down
    // somewhere, with b = k*f the f of where it goes, weighted by what each cell keeps
    auto pcSweepUp = [&](const List<scalarField>& b)
    {
        const label nq = b.size();
        const labelUList& ow = mesh_.owner();
        const labelUList& ne = mesh_.neighbour();
        const cellList& cells = mesh_.cells();
        const label nInt = mesh_.nInternalFaces();
        const scalarField& wI = wPC.primitiveField();
        const scalarField& WI = WPC.primitiveField();
        const scalarField& kI = kPC.primitiveField();
        List<scalarField> X(b);
        PtrList<volScalarField> Xf(nq);
        for (label q = 0; q < nq; ++q)
        {
            Xf.set
            (
                q,
                new volScalarField
                (
                    IOobject("pcX" + Foam::name(q), runTime.timeName(), mesh_),
                    mesh_,
                    dimensionedScalar(dimless, Zero)
                )
            );
            Xf[q].primitiveFieldRef() = b[q];
        }
        scalar tol = 0;
        for (label q = 0; q < nq; ++q) tol = Foam::max(tol, gMax(mag(b[q])));
        tol = 1e-12*tol + VSMALL;
        label sweep = 0;
        scalar change = GREAT;
        for (; sweep < 64 && change > tol; ++sweep)
        {
            List<scalarField> Cq(nq, scalarField(mesh_.nCells(), 0.0));
            List<List<scalarField>> xUsed(nq);
            for (label q = 0; q < nq; ++q)
            {
                xUsed[q].setSize(mesh_.boundary().size());
                forAll(wPC.boundaryField(), patchi)
                {
                    if (!mesh_.boundary()[patchi].coupled()) continue;
                    const scalarField& wP = wPC.boundaryField()[patchi];
                    const labelUList& fc = mesh_.boundary()[patchi].faceCells();
                    xUsed[q][patchi] = Xf[q].boundaryField()[patchi].patchNeighbourField();
                    forAll(fc, i)
                    {
                        if (wP[i] > 0) Cq[q][fc[i]] += wP[i]*xUsed[q][patchi][i];
                    }
                }
            }
            forAllReverse(pcOrder, idx)
            {
                const label c = pcOrder[idx];
                List<scalar> sum(nq, 0.0);
                for (const label facei : cells[c])
                {
                    if (facei >= nInt) continue;
                    const scalar w = wI[facei];
                    label n = -1;
                    if (w > 0 && ow[facei] == c) n = ne[facei];
                    else if (w < 0 && ne[facei] == c) n = ow[facei];
                    if (n < 0) continue;
                    for (label q = 0; q < nq; ++q) sum[q] += mag(w)*X[q][n];
                }
                const scalar pass = (1.0 - kI[c])/WI[c];
                for (label q = 0; q < nq; ++q)
                {
                    X[q][c] = b[q][c] + pass*(sum[q] + Cq[q][c]);
                }
            }
            if (!pcCoupled) break;
            change = 0;
            for (label q = 0; q < nq; ++q)
            {
                Xf[q].primitiveFieldRef() = X[q];
                Xf[q].correctBoundaryConditions();
                forAll(wPC.boundaryField(), patchi)
                {
                    if (!mesh_.boundary()[patchi].coupled()) continue;
                    const scalarField& wP = wPC.boundaryField()[patchi];
                    const scalarField xNb(Xf[q].boundaryField()[patchi].patchNeighbourField());
                    forAll(wP, i)
                    {
                        if (wP[i] > 0)
                        {
                            change = Foam::max(change, mag(xNb[i] - xUsed[q][patchi][i]));
                        }
                    }
                }
            }
            reduce(change, maxOp<scalar>());
        }
        if (pcCoupled && change > tol)
        {
            Info<< "phase change: coupled sweep (up) not converged after " << sweep
                << " sweeps, change " << change << endl;
        }
        return X;
    };


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

    // Raoult-Dalton kinetics: x_i = mole fraction in the liquid (condensed group), y_i = mole fraction in
    // the gas including air; the gas side sees only its partial pressure y_i*p



    // heat capacity TEqn.H applies (updateRhoCv.H, same guard), energy clamp only
    const bool ledgerFloorEOS
    (
        phasedictionary.getOrDefault<bool>("ledgerFloorEOS", false)
    );
    volScalarField rhoCvMix
    (
        IOobject("rhoCvMix", runTime.timeName(), mesh_),
        mesh_,
        dimensionedScalar("rhoCvMix", dimEnergy/dimVolume/dimTemperature, 0.0)
    );
    {
        const dimensionedScalar cZeroCv("cZeroCv", dimDensity, 0.0);
        label i = 0;
        for (const phaseModel& ph : phases_)
        {
            volScalarField cEff(max(c_[i], cZeroCv));
            if (ledgerFloorEOS)
            {
                cEff = max(cEff, ph*ph.thermo().rho());
            }
            rhoCvMix += cEff*ph.thermo().Cv();
            ++i;
        }
    }


    //
    //  molar densities n_i = max(c_i,0)/W_i then mole fractions within each
    //  group. built straight from teh mass-balance

    PtrList<volScalarField> moleFrac(nPhases);
    scalarField nGasI(mesh_.nCells(), 0.0);
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
        nGasI = nGas.primitiveField();

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

    // implicit limit: each cell's gas region (connected cells with beta < 0.5; condensed-side cells join their largest
    // gas neighbour's, two layers deep) and the region's gas [kmol]
    labelList gasRegion(mesh_.nCells(), -1);
    scalarField regN;
    {
        const scalarField& bI = condensate.primitiveField();
        const scalarField& Vc = mesh_.V().field();
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        boolList blocked(mesh_.nFaces(), false);
        forAll(nei, facei)
        {
            blocked[facei] = ((bI[own[facei]] < 0.5) != (bI[nei[facei]] < 0.5));
        }
        forAll(mesh_.boundary(), patchi)
        {
            const fvPatch& fvp = mesh_.boundary()[patchi];
            if (!fvp.coupled()) continue;
            const scalarField bNb(condensate.boundaryField()[patchi].patchNeighbourField());
            const labelUList& fc = fvp.faceCells();
            forAll(fc, i)
            {
                blocked[fvp.start() + i] = ((bI[fc[i]] < 0.5) != (bNb[i] < 0.5));
            }
        }
        const regionSplit rs(mesh_, blocked);
        regN.setSize(rs.nRegions(), 0.0);
        forAll(bI, celli)
        {
            if (bI[celli] < 0.5)
            {
                gasRegion[celli] = rs[celli];
                regN[rs[celli]] += nGasI[celli]*Vc[celli];
            }
        }
        if (Pstream::parRun()) reduce(regN, sumOp<scalarField>());

        volScalarField gR
        (
            IOobject("phaseChangeGasRegion", runTime.timeName(), mesh_),
            mesh_,
            dimensionedScalar(dimless, -1)
        );
        for (label pass = 0; pass < 2; ++pass)
        {
            forAll(gasRegion, celli)
            {
                gR.primitiveFieldRef()[celli] = gasRegion[celli];
            }
            gR.correctBoundaryConditions();
            labelList adopt(gasRegion);
            auto offer = [&](const label c, const label g)
            {
                if (bI[c] < 0.5 || g < 0) return;
                if (adopt[c] < 0 || regN[g] > regN[adopt[c]]) adopt[c] = g;
            };
            forAll(nei, facei)
            {
                offer(own[facei], gasRegion[nei[facei]]);
                offer(nei[facei], gasRegion[own[facei]]);
            }
            forAll(mesh_.boundary(), patchi)
            {
                if (!mesh_.boundary()[patchi].coupled()) continue;
                const scalarField gNb(gR.boundaryField()[patchi].patchNeighbourField());
                const labelUList& fc = mesh_.boundary()[patchi].faceCells();
                forAll(fc, i)
                {
                    offer(fc[i], label(std::round(gNb[i])));
                }
            }
            gasRegion = adopt;
        }
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


                        // density floor for dVspec = 1/rhoVap in PCR: a cell floored at 1e-3 gave the pressure equation a
                        // coefficient 100x a healthy cell's; 0.01 caps 1/rho at 100 (phaseChangeRhoFloor)
                        const dimensionedScalar rhoFloorPC
                        (
                            "rhoFloorPC",
                            dimDensity,
                            phasedictionary.getOrDefault<scalar>
                            (
                                "phaseChangeRhoFloor", 0.01
                            )
                        );
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

                        // the two sides of the interface: liquid composition from the
                        // condensed cells, vapour partial pressure from the gas cells
                        const volScalarField xLiqI(interfaceState(xLiq, condensate));
                        // vapour partial pressure the surface works against: that of the gas its vapour goes to (or comes from),
                        // weighted as the sweep lays it down; a cell whose vapour goes nowhere takes its neighbours' face mean
                        volScalarField pVapI
                        (
                            IOobject("pVapI." + liq.name(), runTime.timeName(), mesh_),
                            interfaceState(yVap*p_, scalar(1) - condensate)
                        );
                        scalarField X0;
                        {
                            List<scalarField> b(2);
                            b[0] = kPC.primitiveField();
                            b[1] = kPC.primitiveField()*(yVap*p_)().primitiveField();
                            const List<scalarField> X(pcSweepUp(b));
                            X0 = X[0];
                            const scalarField& WI = WPC.primitiveField();
                            scalarField& pI = pVapI.primitiveFieldRef();
                            forAll(pI, celli)
                            {
                                if (WI[celli] > 0 && X[0][celli] > VSMALL)
                                {
                                    pI[celli] = X[1][celli]/X[0][celli];
                                }
                            }
                            pVapI.correctBoundaryConditions();
                        }

                        // Hertz-Knudsen flux [kg/m2/s], positive evaporates: x*Psat against the vapour partial
                        // pressure there, no Schrage drift factor (review: Persad & Ward 2016,
                        // https://doi.org/10.1021/acs.chemrev.5b00511)
                        const volScalarField jHK
                        (
                            "jHK." + liq.name(),
                            Foam::sqrt
                            (
                                (Wliq/1000.0)/(2.0*M_PI*gasconstant*Tsafe)
                            )
                           *(xLiqI*Psat - pVapI)
                        );

                        // the anchor rate for this pair [kg/m3/s], positive evaporates:
                        // the flux times the accommodation and the interface area density
                        volScalarField r
                        (
                            IOobject("r" + liq.name(), runTime.timeName(), mesh_),
                            jHK*(pos0(jHK)*sigmaEvap + neg(jHK)*sigmaCond)*aI
                        );

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

                        // implicit limit (Zenz et al. 2024, Eq. 67, https://doi.org/10.1016/j.compfluid.2023.106109, with the vapour's volume): s = this step's rise
                        // in the receivers' y*p (own gas for y, gas region at fixed volume for p) over x*Psat - pVapI; r /= s where s > 1, so no step overshoots
                        scalar minFac = 1;
                        label nImp = 0;
                        {
                            const scalarField& Vc = mesh_.V().field();
                            const scalarField& kk = kPC.primitiveField();
                            const scalarField& yI = yVap.primitiveField();
                            const scalarField& pI = p_.primitiveField();
                            const volScalarField Wvap(vap.thermo().W());
                            const scalarField& WvI = Wvap.primitiveField();

                            const scalarField rV0(r.primitiveField()*Vc);
                            const List<scalarField> Yd(pcSweepDown(List<scalarField>(1, rV0)));
                            const scalarField dm(kk*Yd[0]*dt);

                            volScalarField dpV
                            (
                                IOobject("dpV." + liq.name(), runTime.timeName(), mesh_),
                                mesh_,
                                dimensionedScalar(dimPressure, Zero)
                            );
                            scalarField& dI = dpV.primitiveFieldRef();
                            forAll(dI, j)
                            {
                                if (dm[j] == 0) continue;
                                const scalar Nown = Foam::max(nGasI[j]*Vc[j], VSMALL);
                                const scalar Nreg =
                                    (gasRegion[j] >= 0) ? Foam::max(regN[gasRegion[j]], Nown) : Nown;
                                dI[j] = dm[j]/WvI[j]*pI[j]*((1 - yI[j])/Nown + yI[j]/Nreg);
                            }
                            dpV.correctBoundaryConditions();
                            const scalarField kdp(kk*dI);
                            const List<scalarField> Xd(pcSweepUp(List<scalarField>(1, kdp)));

                            const scalarField& xI = xLiqI.primitiveField();
                            const scalarField& PsI = Psat.primitiveField();
                            const scalarField& pVI = pVapI.primitiveField();
                            scalarField& rI = r.primitiveFieldRef();
                            forAll(rI, a)
                            {
                                const scalar drive = xI[a]*PsI[a] - pVI[a];
                                if (rI[a] == 0 || drive == 0 || !(X0[a] > VSMALL)) continue;
                                const scalar s = Xd[0][a]/(X0[a]*drive);
                                if (s <= 1) continue;
                                rI[a] /= s;
                                minFac = Foam::min(minFac, 1/s);
                                ++nImp;
                            }
                        }
                        reduce(minFac, minOp<scalar>());
                        reduce(nImp, sumOp<label>());

                        
                        // saturation temperature consistent with the composition, x*Psat(T) = y*p (zero driving force), so the
                        // energy clamp below brakes to the kinetics' own fixed point: Raoult with the integrated Clausius-
                        // Clapeyron law inverted, Tsat = Tboil/(1 - ln(y*p/(x*P0))/K), K = M*Lvap/(R*Tboil)
                        
                        const volScalarField Kfield //for above tect relation///
                        (
                            ((Wliq/1000.0)*Lvap)/(Tboil*gasconstant)
                        );
                        const volScalarField satArg
                        (
                            max
                            (
                                pVapI
                               /(
                                    max
                                    (
                                        xLiqI,
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



                        label nClampE = 0, nClampM = 0, nClampV = 0;

                        // anchor limits: the latent heat is exchanged in the anchor, so it
                        // may not drive the anchor past Tsat this step, and evaporation
                        // takes at most half the condensed mass there
                        {
                            scalarField& rI = r.primitiveFieldRef();
                            const scalarField& TI      = Tsafe.primitiveField();
                            const scalarField& TsatI   = Tsat.primitiveField();
                            const scalarField& rhoCvMixI = rhoCvMix.primitiveField();
                            // clamp budget: the anchor's whole cooling per kg is Lvap
                            const scalar       LvapV   = Lvap.value();
                            const scalarField& cLiqI = c_[li].primitiveField();
                            forAll(rI, celli)//LIMITING SO CODE DOESNT BLOW UP
                            {
                                scalar& rc = rI[celli];
                                const scalar rcKinetic = rc;

                                //  energy limit: dont let latent heat drive
                                //  T past Tsat this step (eFac = rhoCv/(dt*Lvap))
                                const scalar eFac =
                                    Foam::max(rhoCvMixI[celli], 0.0)*rDeltaT/LvapV;
                                const scalar headUp =
                                    Foam::max(TsatI[celli] - TI[celli], 0.0);
                                const scalar headDn =
                                    Foam::max(TI[celli] - TsatI[celli], 0.0);
                                const scalar rMinE = -headUp*eFac;  // cond cap (<=0)
                                const scalar rMaxE =  headDn*eFac;  // evap cap (>=0)
                                // count only clamps that matter: a trace rate clamped to a smaller trace is noise (14k cells at 1e-51)
                                if (rc > rMaxE) { rc = rMaxE; if (mag(rcKinetic) > 1e-3) ++nClampE; }
                                if (rc < rMinE) { rc = rMinE; if (mag(rcKinetic) > 1e-3) ++nClampE; }

                                //  mass limit: half of whats there per step,
                                //  max. NB the positivity limiter budget below
                                //  RELIES on this 0.5, dont raise it
                                const scalar rmaxEvap =
                                    0.5*Foam::max(cLiqI[celli], 0.0)*rDeltaT;
                                if (rc >  rmaxEvap) { rc =  rmaxEvap; ++nClampM; }
                            }
                        }

                        const scalarField& V = mesh_.V().field();
                        const volScalarField& TT = T_;
                        const dimensionedScalar qZero("qZero", dimMass/dimTime/dimArea, Zero);

                        // per unit face share [kg/s/m2]: q = rate*V/W leaves through each unit
                        volScalarField qPC
                        (
                            IOobject("qPC." + liq.name(), runTime.timeName(), mesh_),
                            mesh_,
                            qZero
                        );
                        auto setQ = [&](const scalarField& rate)
                        {
                            const scalarField& WI = WPC.primitiveField();
                            scalarField& qI = qPC.primitiveFieldRef();
                            forAll(qI, celli)
                            {
                                qI[celli] =
                                    (WI[celli] > 0) ? rate[celli]*V[celli]/WI[celli] : 0.0;
                            }
                            qPC.correctBoundaryConditions();
                        };

                        // the evaporating liquid goes from the top: an anchor's share moves into the neighbour with less condensed
                        // phase, as far as it has liquid; the heat still comes from the anchor and the vapour is born at its T.
                        // Taking liquid from under a part-filled cell left a gas pocket that boiled the film apart (test 07)
                        setQ(r.primitiveField());
                        const volScalarField qEv("qEv", max(qPC, qZero));
                        const scalarField offer(pcReceive(qEv));
                        volScalarField fShift
                        (
                            IOobject("fShift." + liq.name(), runTime.timeName(), mesh_),
                            mesh_,
                            dimensionedScalar("fShift", dimless, Zero)
                        );
                        {
                            const scalarField& rI = r.primitiveField();
                            const scalarField& WI = WPC.primitiveField();
                            const scalarField& cLiqI = c_[li].primitiveField();
                            scalarField& fI = fShift.primitiveFieldRef();
                            forAll(fI, celli)
                            {
                                if (!(WI[celli] > 0) || !(offer[celli] > 0)) continue;
                                const scalar budget =
                                    (0.5*Foam::max(cLiqI[celli], 0.0)*rDeltaT
                                  - Foam::max(rI[celli], 0.0))*V[celli];
                                fI[celli] =
                                    Foam::min(Foam::max(budget, 0.0)/offer[celli], 1.0);
                            }
                        }
                        fShift.correctBoundaryConditions();

                        const scalarField& qEvI = qEv.primitiveField();
                        const scalarField& fI = fShift.primitiveField();
                        const scalarField& TI = TT.primitiveField();
                        const scalarField give(qEvI*pcSendWeighted(fShift));
                        const scalarField take(fI*offer);
                        const scalarField takeT(fI*pcReceive(volScalarField(qEv*TT)));

                        // liquid made vapour in each cell [kg/m3/s], and the temperature that
                        // vapour is born at (its anchors' temperature)
                        scalarField rS(r.primitiveField() + (take - give)/V);
                        volScalarField TbS("TbS." + liq.name(), TT);
                        {
                            const scalarField& rI = r.primitiveField();
                            scalarField& TbI = TbS.primitiveFieldRef();
                            forAll(TbI, celli)
                            {
                                const scalar own = Foam::max(rI[celli]*V[celli] - give[celli], 0.0);
                                if (own + take[celli] > 0)
                                {
                                    TbI[celli] =
                                        (own*TI[celli] + takeT[celli])/(own + take[celli]);
                                }
                            }
                        }
                        TbS.correctBoundaryConditions();

                        // vapour this step [kg/s] per cell, own (rS, born at TbS) and passed down the shares: evaporated vapour
                        // carries its birth velocity and T, condensing vapour the giving cell's. 0 evaporated, 1-3 its momentum,
                        // 4 mass times birth T, 5 condensing (<= 0)
                        const vectorField& UI = U_.primitiveField();
                        const scalarField& kI = kPC.primitiveField();
                        List<scalarField> Yv;
                        {
                            List<scalarField> own(6, scalarField(mesh_.nCells(), 0.0));
                            const scalarField& TbI = TbS.primitiveField();
                            forAll(rS, celli)
                            {
                                const scalar m = rS[celli]*V[celli];
                                if (m > 0)
                                {
                                    own[0][celli] = m;
                                    own[1][celli] = m*UI[celli].x();
                                    own[2][celli] = m*UI[celli].y();
                                    own[3][celli] = m*UI[celli].z();
                                    own[4][celli] = m*TbI[celli];
                                }
                                else
                                {
                                    own[5][celli] = m;
                                }
                            }
                            Yv = pcSweepDown(own);
                        }

                        // limits: condensation takes at most half a cell's vapour per step, and the volume made or removed stays
                        // under etaVol of the cell. sPCm (mass limit alone) is what the dt control resolves, so the cap never
                        // truncates evaporative cooling. What a limit holds back is refused at the cells it came from
                        volScalarField sPC
                        (
                            IOobject("sPC." + liq.name(), runTime.timeName(), mesh_),
                            mesh_,
                            dimensionedScalar("sPC", dimless, 1.0)
                        );
                        volScalarField sPCm("sPCm." + liq.name(), sPC);
                        {
                            const scalarField& cVapI = c_[vi].primitiveField();
                            const scalarField& rhoLiqI = rhoLiq.primitiveField();
                            const scalarField& rhoVapI = rhoVap.primitiveField();
                            scalarField& sI = sPC.primitiveFieldRef();
                            scalarField& smI = sPCm.primitiveFieldRef();
                            forAll(sI, celli)
                            {
                                const scalar R = kI[celli]*(Yv[0][celli] + Yv[5][celli]);
                                if (R == 0) continue;

                                scalar sm = 1.0;
                                const scalar capM =
                                    0.5*Foam::max(cVapI[celli], 0.0)*V[celli]*rDeltaT;
                                if (-R > capM) { sm = capM/(-R); ++nClampM; }

                                // volume limit etaVol of the cell per step, a rate cap ~1/dt: the dt control's maxPhaseChangeFraction
                                // must stay below it (see the etaVol read)
                                const scalar dVspec =
                                    Foam::max
                                    (
                                        1.0/rhoVapI[celli] - 1.0/rhoLiqI[celli],
                                        VSMALL
                                    );
                                const scalar capV = etaVol*V[celli]*rDeltaT/dVspec;
                                scalar sv = 1.0;
                                if (mag(R)*sm > capV) { sv = capV/(mag(R)*sm); ++nClampV; }

                                smI[celli] = sm;
                                sI[celli] = sm*sv;
                            }
                        }
                        sPC.correctBoundaryConditions();
                        sPCm.correctBoundaryConditions();

                        // up the shares: the part of each cell's vapour that is laid down (0,
                        // with both limits; 1, the mass limit alone), and for condensation
                        // the momentum (2-4) and temperature (5) of the vapour it draws
                        List<scalarField> Av;
                        {
                            const scalarField& sI = sPC.primitiveField();
                            const scalarField& smI = sPCm.primitiveField();
                            List<scalarField> b(6);
                            b[0] = sI*kI;
                            b[1] = smI*kI;
                            b[2] = sI*kI*UI.component(vector::X);
                            b[3] = sI*kI*UI.component(vector::Y);
                            b[4] = sI*kI*UI.component(vector::Z);
                            b[5] = sI*kI*TI;
                            Av = pcSweepUp(b);
                        }
                        volScalarField accPC
                        (
                            IOobject("accPC." + liq.name(), runTime.timeName(), mesh_),
                            mesh_,
                            dimensionedScalar("accPC", dimless, Zero)
                        );
                        accPC.primitiveFieldRef() = Av[0];
                        accPC.correctBoundaryConditions();

                        // paired rates [kg/m3/s]: liquid made vapour in its cell (rL), vapour laid down (rV) and the anchor's
                        // share (rH) that pays the latent heat; the totals pair exactly
                        volScalarField rL
                        (
                            IOobject("rL." + liq.name(), runTime.timeName(), mesh_),
                            mesh_,
                            dimensionedScalar("rL", dimDensity/dimTime, Zero)
                        );
                        volScalarField rV("rV." + liq.name(), rL);
                        volScalarField rH("rH." + liq.name(), rL);
                        volScalarField rLpre("rLpre." + liq.name(), rL);
                        volScalarField rVpre("rVpre." + liq.name(), rL);
                        {
                            const scalarField& rI = r.primitiveField();
                            const scalarField& sI = sPC.primitiveField();
                            const scalarField& smI = sPCm.primitiveField();
                            const scalarField passedOn
                            (
                                qEvI*pcSendWeighted(volScalarField(fShift*accPC))
                            );
                            scalarField& rLI = rL.primitiveFieldRef();
                            scalarField& rVI = rV.primitiveFieldRef();
                            scalarField& rHI = rH.primitiveFieldRef();
                            scalarField& rLpI = rLpre.primitiveFieldRef();
                            scalarField& rVpI = rVpre.primitiveFieldRef();
                            forAll(rLI, celli)
                            {
                                const scalar kept = kI[celli]*(Yv[0][celli] + Yv[5][celli]);
                                rLI[celli] = rS[celli]*Av[0][celli];
                                rLpI[celli] = rS[celli]*Av[1][celli];
                                rVI[celli] = sI[celli]*kept/V[celli];
                                rVpI[celli] = smI[celli]*kept/V[celli];
                                rHI[celli] =
                                    ((rI[celli]*V[celli] - give[celli])*Av[0][celli]
                                  + passedOn[celli])/V[celli];
                            }
                        }

                        SuMass[li] -= rL;    // condensed side loses (gains on condensation)
                        SuMass[vi] += rV;    // gas side gains, the same mass

                        // phase-change dilatation as applied (contraction where the liquid
                        // goes, expansion where the vapour arrives), and before the volume cap
                        // for the dt control
                        PCR += rV/rhoVap - rL/rhoLiq;
                        pcGasDil_ += rV/rhoVap;
                        pcCondDil_ -= rL/rhoLiq;
                        if (vDotDtPtr)
                        {
                            *vDotDtPtr += rVpre/rhoVap - rLpre/rhoLiq;
                        }

                        // heat via -mass_dot in TEqn [W/m3]: the anchor gives the enthalpy of vaporisation (LvapE + p*dv = Lvap)
                        // and heats the liquid it evaporates in a neighbour; each cell gets back the p*dV its work term takes for
                        // the phase-change volume; vapour is laid down at its birth T (else a receiver overheats, test 09)
                        {
                            const volScalarField CvV(vap.thermo().Cv());
                            const volScalarField CvL(liq.thermo().Cv());
                            const scalarField shiftW(pcSendWeighted(volScalarField(fShift*accPC)));
                            const scalarField shiftWT
                            (
                                pcSendWeighted(volScalarField(fShift*accPC*TT))
                            );
                            const scalarField& sI = sPC.primitiveField();
                            const scalarField& CvVI = CvV.primitiveField();
                            const scalarField& CvLI = CvL.primitiveField();
                            const scalarField& pI = p_.primitiveField();
                            const scalarField& rhoLiqI = rhoLiq.primitiveField();
                            const scalarField& rhoVapI = rhoVap.primitiveField();
                            const scalarField& rLI = rL.primitiveField();
                            const scalarField& rVI = rV.primitiveField();
                            const scalarField& rHI = rH.primitiveField();
                            const scalar LvapV = Lvap.value();
                            scalarField& mdI = massdotterm->primitiveFieldRef();
                            forAll(mdI, celli)
                            {
                                const scalar mCond = Foam::min(rS[celli], 0.0)*V[celli];
                                const scalar heatIn =
                                    CvVI[celli]*sI[celli]*kI[celli]
                                   *(Yv[4][celli] - TI[celli]*Yv[0][celli])
                                  - CvVI[celli]*mCond*(Av[5][celli] - TI[celli]*Av[0][celli])
                                  - CvLI[celli]*qEvI[celli]*(TI[celli]*shiftW[celli] - shiftWT[celli]);
                                mdI[celli] +=
                                    rHI[celli]*LvapV
                                  + rLI[celli]*pI[celli]/rhoLiqI[celli]
                                  - rVI[celli]*pI[celli]/rhoVapI[celli]
                                  - heatIn/V[celli];
                            }
                        }

                        // momentum: the mass carries its velocity across (vapour is born with its liquid's, condensate joins
                        // with the vapour's), so a moving melt hands its velocity to its vapour
                        {
                            const scalarField& sI = sPC.primitiveField();
                            vectorField& mI = pcMomentum_.primitiveFieldRef();
                            forAll(mI, celli)
                            {
                                const scalar Vc = V[celli];
                                const scalar f = sI[celli]*kI[celli];
                                // evaporation: out of the liquid's cell, into where it is laid
                                const scalar mEv = Foam::max(rS[celli], 0.0)*Vc*Av[0][celli];
                                mI[celli] +=
                                    (
                                        f*vector(Yv[1][celli], Yv[2][celli], Yv[3][celli])
                                      - mEv*UI[celli]
                                    )/Vc;
                                // condensation: out of the cells giving the vapour, into the
                                // condensing one
                                const scalar mCond = Foam::min(rS[celli], 0.0)*Vc;
                                mI[celli] +=
                                    (
                                        f*Yv[5][celli]*UI[celli]
                                      - mCond*vector(Av[2][celli], Av[3][celli], Av[4][celli])
                                    )/Vc;
                            }
                        }

                        // (no implicit pressure coupling kP: it added pEqn volume no mass
                        // backed, and the work term turned it into energy; 1D test 08)

                        reduce(nClampE, sumOp<label>());
                        reduce(nClampM, sumOp<label>());
                        reduce(nClampV, sumOp<label>());

                        const scalarField& rBud = rL.primitiveField();
                        const scalar evapKgPerS =
                            gSum(pos(rBud)*rBud*mesh_.V().field());
                        const scalar condKgPerS =
                            gSum(neg(rBud)*rBud*mesh_.V().field());
                        const scalar pairErr =
                            gSum((rBud - rV.primitiveField())*mesh_.V().field());

                        Info<< "    (" << liq.name() << "," << vap.name()
                            << ")  Tboil=" << Tboil.value()
                            << "  max rate magnitude=" << max(mag(rL)).value()
                            << " kg/m3/s"
                            << "  clamped cells ENERGY/MASS/VOLUME to maintain boundedness..... probably decrease timestep or pray to the numerical gods: "
                            << nClampE << "/" << nClampM << "/" << nClampV
                            << endl;
                        Info<< "        integrated evap = " << evapKgPerS
                            << " kg/s, cond = " << condKgPerS
                            << " kg/s (cond is <= 0), pairing " << pairErr << " kg/s"
                            << endl;
                        Info<< "        implicit limit: " << nImp << " cell(s), smallest factor " << minFac << endl;
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

    // phiPC [m3/s]: each cell's phase-change volume plus what reaches it from its condensed side, passed on where
    // beta drops. The condensed group rides phi - phiPC, as in https://doi.org/10.1016/j.jcp.2019.108954
    // (Palmore & Desjardins 2019) and https://doi.org/10.1016/j.jcp.2020.109251 (Scapin et al. 2020)
    surfaceScalarField phiPC
    (
        IOobject("phiPC", runTime.timeName(), mesh_),
        mesh_,
        dimensionedScalar("phiPC", dimVolume/dimTime, Zero)
    );
    {
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();
        const cellList& cells = mesh_.cells();
        const label nIntFaces = mesh_.nInternalFaces();
        const scalarField& wI = wPC.primitiveField();
        const scalarField& WI = WPC.primitiveField();
        const scalarField& bI = condensate.primitiveField();
        const scalarField Q(PCR.primitiveField()*mesh_.V().field());

        // outflow per unit share [m/s]
        volScalarField gPC
        (
            IOobject("gPC", runTime.timeName(), mesh_),
            mesh_,
            dimensionedScalar("gPC", dimVelocity, Zero)
        );
        scalarField& gI = gPC.primitiveFieldRef();

        // the cells that pass volume on, by decreasing beta: everything upstream of a
        // cell is final when the cell is visited, so one sweep is exact on one processor
        labelList order;
        {
            DynamicList<label> band;
            forAll(WI, celli)
            {
                if (WI[celli] > 0) band.append(celli);
            }
            scalarField key(band.size());
            forAll(band, i) key[i] = -bI[band[i]];
            labelList idx;
            sortedOrder(key, idx);
            order.setSize(band.size());
            forAll(idx, i) order[i] = band[idx[i]];
        }

        bool anyCoupled = false;
        forAll(mesh_.boundary(), patchi)
        {
            if (mesh_.boundary()[patchi].coupled()) anyCoupled = true;
        }
        reduce(anyCoupled, orOp<bool>());

        // across coupled faces the inflow from the other side is iterated until it no
        // longer changes (a chain crossing k processor faces needs k + 1 sweeps)
        const scalar tolPC = 1e-12*gSum(mag(Q)) + VSMALL;
        label sweep = 0;
        scalar change = GREAT;
        for (; sweep < 64 && change > tolPC; ++sweep)
        {
            scalarField I(mesh_.nCells(), 0.0);
            List<scalarField> gNbUsed(mesh_.boundary().size());
            forAll(wPC.boundaryField(), patchi)
            {
                if (!mesh_.boundary()[patchi].coupled()) continue;
                const scalarField& wP = wPC.boundaryField()[patchi];
                const labelUList& fc = mesh_.boundary()[patchi].faceCells();
                gNbUsed[patchi] = gPC.boundaryField()[patchi].patchNeighbourField();
                const scalarField& gNb = gNbUsed[patchi];
                forAll(fc, i)
                {
                    if (wP[i] < 0) I[fc[i]] -= gNb[i]*wP[i];
                }
            }

            for (const label celli : order)
            {
                gI[celli] = (Q[celli] + I[celli])/WI[celli];
                for (const label facei : cells[celli])
                {
                    if (facei >= nIntFaces) continue;
                    const scalar w = wI[facei];
                    if (w > 0 && own[facei] == celli) I[nei[facei]] += gI[celli]*w;
                    else if (w < 0 && nei[facei] == celli) I[own[facei]] -= gI[celli]*w;
                }
            }

            if (!anyCoupled) break;

            gPC.correctBoundaryConditions();
            change = 0;
            forAll(wPC.boundaryField(), patchi)
            {
                if (!mesh_.boundary()[patchi].coupled()) continue;
                const scalarField& wP = wPC.boundaryField()[patchi];
                const scalarField gNb(gPC.boundaryField()[patchi].patchNeighbourField());
                forAll(wP, i)
                {
                    change = Foam::max(change, mag((gNb[i] - gNbUsed[patchi][i])*wP[i]));
                }
            }
            reduce(change, maxOp<scalar>());
        }
        if (anyCoupled && change > tolPC)
        {
            Info<< "phiPC: coupled sweep not converged after " << sweep
                << " sweeps, change " << change << " m3/s" << endl;
        }

        forAll(own, facei)
        {
            const scalar w = wI[facei];
            if (w > 0) phiPC[facei] = gI[own[facei]]*w;
            else if (w < 0) phiPC[facei] = gI[nei[facei]]*w;
        }
        forAll(phiPC.boundaryField(), patchi)
        {
            scalarField& pP = phiPC.boundaryFieldRef()[patchi];
            pP = 0.0;
            if (!mesh_.boundary()[patchi].coupled()) continue;
            const scalarField& wP = wPC.boundaryField()[patchi];
            const labelUList& fc = mesh_.boundary()[patchi].faceCells();
            const scalarField gNb(gPC.boundaryField()[patchi].patchNeighbourField());
            forAll(fc, i)
            {
                if (wP[i] > 0) pP[i] = gI[fc[i]]*wP[i];
                else if (wP[i] < 0) pP[i] = gNb[i]*wP[i];
            }
        }
    }

    surfaceScalarField phiL("phiL", phi_ - phiPC);

    // raw un-limited beta flux
    surfaceScalarField betaPhiCorr
    (
        "betaPhiCorr",
        fvc::flux(phiL, beta, alphaScheme)
    );

    // interFoam style compresion flux normal to the interface, keeps it sharp
    {
        surfaceScalarField phic(mag(phiL)/mesh_.magSf());
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
        phiL,
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
        betaPhiCorr + upwind<scalar>(mesh_, phiL).flux(beta)
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
    rhoCvPhi_ = dimensionedScalar(dimEnergy/dimTemperature/dimTime, Zero);
    rhoCvTPhi_ = dimensionedScalar(dimEnergy/dimTime, Zero);

    // betaPhi + gammaPhi == phi_: the volume the ledger moves this step, TEqn's p*dV
    divPhiLedger_ = fvc::div(phi_);


    //  Y REFORMULATION of the species transport (attempt no.3, this one works)



    // receiver packing cap, off unless cPackTol > 0 (legacy 1.02): it throttles
    // the mass of a face but not its volume, so heavy or cold inflow arrives as volume without mass
    const scalar cPackTol =
        phasedictionary.getOrDefault<scalar>("cPackTol", 0.0);

    // --- group totals, group sources, and the packing cap density ---------
         //  cap = densest species of the group at local (p,T), tightened below to
        //  the compositon weighted mixture density where theres actual content
    volScalarField cCondTot
    (
        IOobject("cCondTot", runTime.timeName(), mesh_),
        mesh_, dimensionedScalar("z", dimDensity, 0.0)
    );
    volScalarField cGasTot
    (
        IOobject("cGasTot", runTime.timeName(), mesh_),
        mesh_, dimensionedScalar("z", dimDensity, 0.0)
    );
    volScalarField SuCondSum
    (
        IOobject("SuCondSum", runTime.timeName(), mesh_),
        mesh_, dimensionedScalar("z", dimensionSet(1, -3, -1, 0, 0), 0.0)
    );
    volScalarField SuGasSum
    (
        IOobject("SuGasSum", runTime.timeName(), mesh_),
        mesh_, dimensionedScalar("z", dimensionSet(1, -3, -1, 0, 0), 0.0)
    );
    volScalarField rhoCondCap
    (
        IOobject("rhoCondCap", runTime.timeName(), mesh_),
        mesh_, dimensionedScalar("z", dimDensity, 0.0)
    );
    volScalarField rhoGasCap
    (
        IOobject("rhoGasCap", runTime.timeName(), mesh_),
        mesh_, dimensionedScalar("z", dimDensity, 0.0)
    );

    volScalarField vCondImp
    (
        IOobject("vCondImp", runTime.timeName(), mesh_),
        mesh_, dimensionedScalar("z", dimless, 0.0)
    );
    volScalarField vGasImp
    (
        IOobject("vGasImp", runTime.timeName(), mesh_),
        mesh_, dimensionedScalar("z", dimless, 0.0)
    );

    {
        const dimensionedScalar rhoFloorY("rhoFloorY", dimDensity, 1e-3);
        label i = 0;
        for (phaseModel& ph : phases_)
        {
            const volScalarField cPos(max(c_[i], dimensionedScalar("z", dimDensity, 0.0)));
            if (isGas[i])
            {
                cGasTot   += cPos;
                vGasImp   += cPos/max(ph.thermo().rho(), rhoFloorY);
                SuGasSum  += SuMass[i];
                rhoGasCap  = max(rhoGasCap, ph.thermo().rho());
            }
            else
            {
                cCondTot  += cPos;
                vCondImp  += cPos/max(ph.thermo().rho(), rhoFloorY);
                SuCondSum += SuMass[i];
                rhoCondCap = max(rhoCondCap, ph.thermo().rho());
            }
            ++i;
        }
    }

    // tighten the packing caaps 
    {
        const scalar vTinyCap = 1e-12;

        scalarField& capCI = rhoCondCap.primitiveFieldRef();
        scalarField& capGI = rhoGasCap.primitiveFieldRef();
        const scalarField& cTC = cCondTot.primitiveField();
        const scalarField& cTG = cGasTot.primitiveField();
        const scalarField& vIC = vCondImp.primitiveField();
        const scalarField& vIG = vGasImp.primitiveField();

        forAll(capCI, celli)
        {
            if (vIC[celli] > vTinyCap)
            {
                capCI[celli] =
                    Foam::min(capCI[celli], cTC[celli]/vIC[celli]);
            }
            if (vIG[celli] > vTinyCap)
            {
                capGI[celli] =
                    Foam::min(capGI[celli], cTG[celli]/vIG[celli]);
            }
        }

        rhoCondCap.correctBoundaryConditions();
        rhoGasCap.correctBoundaryConditions();
    }

 
    // group density of the donor's own content, uncapped: an over-full donor sheds its excess
    const volScalarField psiCond
    (
        "psiCond",
        max
        (
            cCondTot/max(condensate, betaSmall),
            dimensionedScalar("zeroRho", dimDensity, 0.0)
        )
    );
    const volScalarField psiGas
    (
        "psiGas",
        max
        (
            cGasTot/max(scalar(1) - condensate, betaSmall),
            dimensionedScalar("zeroRho", dimDensity, 0.0)
        )
    );

    surfaceScalarField FCond("FCond", upwind<scalar>(mesh_, betaPhi).flux(psiCond));
    surfaceScalarField FGas ("FGas",  upwind<scalar>(mesh_, gammaPhi).flux(psiGas));

    // physical OUTFLOW faces: upwind takes the patch value (c_b = alpha_b*rho_b, the inflow composition); the donor is the cell
    auto donorOutflow = [&]
    (
        surfaceScalarField& F,
        const surfaceScalarField& dirFlux,
        const volScalarField& psi
    )
    {
        auto& Fbf = F.boundaryFieldRef();
        forAll(Fbf, patchi)
        {
            if (Fbf[patchi].coupled()) continue;

            fvsPatchScalarField& pF = Fbf[patchi];
            const fvsPatchScalarField& pDir = dirFlux.boundaryField()[patchi];
            const scalarField psiIn(psi.boundaryField()[patchi].patchInternalField());
            forAll(pF, bfacei)
            {
                if (pDir[bfacei] > 0)
                {
                    pF[bfacei] = pDir[bfacei]*psiIn[bfacei];
                }
            }
        }
    };
    donorOutflow(FCond, betaPhi, psiCond);
    donorOutflow(FGas, gammaPhi, psiGas);

    // face value from the donor of F: upwind inside and on coupled patches, the cell on physical outflow faces
    auto donorValue = [&]
    (
        const surfaceScalarField& F,
        const volScalarField& vf
    )
    {
        tmp<surfaceScalarField> tvf(upwind<scalar>(mesh_, F).interpolate(vf));
        auto& vfBf = tvf.ref().boundaryFieldRef();
        forAll(vfBf, patchi)
        {
            if (vfBf[patchi].coupled()) continue;

            const fvsPatchScalarField& pF = F.boundaryField()[patchi];
            const scalarField vIn(vf.boundaryField()[patchi].patchInternalField());
            forAll(vfBf[patchi], bfacei)
            {
                if (pF[bfacei] > 0)
                {
                    vfBf[patchi][bfacei] = vIn[bfacei];
                }
            }
        }
        return tvf;
    };


    volScalarField lamOutCond
    (
        IOobject("lamOutCond", runTime.timeName(), mesh_),
        mesh_, dimensionedScalar("one", dimless, 1.0)
    );
    volScalarField lamInCond
    (
        IOobject("lamInCond", runTime.timeName(), mesh_),
        mesh_, dimensionedScalar("one", dimless, 1.0)
    );
    volScalarField lamOutGas
    (
        IOobject("lamOutGas", runTime.timeName(), mesh_),
        mesh_, dimensionedScalar("one", dimless, 1.0)
    );
    volScalarField lamInGas
    (
        IOobject("lamInGas", runTime.timeName(), mesh_),
        mesh_, dimensionedScalar("one", dimless, 1.0)
    );

    for (int grp = 0; grp < 2; ++grp)
    {
        surfaceScalarField& Fg   = (grp == 0) ? FCond      : FGas;
        volScalarField& lamO     = (grp == 0) ? lamOutCond : lamOutGas;
        volScalarField& lamN     = (grp == 0) ? lamInCond  : lamInGas;
        const scalarField& cTotI = ((grp == 0) ? cCondTot  : cGasTot).primitiveField();
        const scalarField& SuI   = ((grp == 0) ? SuCondSum : SuGasSum).primitiveField();
        const scalarField& capI  = ((grp == 0) ? rhoCondCap : rhoGasCap).primitiveField();

        const scalarField& V = mesh_.V().field();
        const labelUList& own = mesh_.owner();
        const labelUList& nei = mesh_.neighbour();

        scalarField outflow(mesh_.nCells(), 0.0);
        scalarField inflow(mesh_.nCells(), 0.0);

        {
            const scalarField& FgI = Fg.primitiveField();
            forAll(FgI, facei)
            {
                if (FgI[facei] > 0)
                {
                    outflow[own[facei]] += FgI[facei];
                    inflow[nei[facei]]  += FgI[facei];
                }
                else
                {
                    outflow[nei[facei]] -= FgI[facei];
                    inflow[own[facei]]  -= FgI[facei];
                }
            }
        }

        forAll(Fg.boundaryField(), patchi)
        {
            const fvsPatchScalarField& pF = Fg.boundaryField()[patchi];
            const labelUList& fc = pF.patch().faceCells();
            forAll(pF, bfacei)
            {
                if (pF[bfacei] > 0)
                {
                    outflow[fc[bfacei]] += pF[bfacei];
                }
                else
                {
                    inflow[fc[bfacei]] -= pF[bfacei];
                }
            }
        }

        scalarField& lamOI = lamO.primitiveFieldRef();
        forAll(lamOI, celli)
        {
            if (outflow[celli] > VSMALL)
            {
                const scalar budget =
                    Foam::max(cTotI[celli], scalar(0))*V[celli]*rDeltaT
                  + Foam::min(SuI[celli], scalar(0))*V[celli];

                lamOI[celli] = Foam::min
                (
                    scalar(1),
                    Foam::max(budget, scalar(0))/outflow[celli]
                );
            }
        }

        scalarField& lamNI = lamN.primitiveFieldRef();
        forAll(lamNI, celli)
        {
            if (cPackTol > 0 && inflow[celli] > VSMALL)
            {
                const scalar headroom =
                    (
                        cPackTol*capI[celli]
                      - Foam::max(cTotI[celli], scalar(0))
                    )*V[celli]*rDeltaT
                  - Foam::max(SuI[celli], scalar(0))*V[celli]
                  + lamOI[celli]*outflow[celli];

                lamNI[celli] = Foam::min
                (
                    scalar(1),
                    Foam::max(headroom, scalar(0))/inflow[celli]
                );
            }
        }

        lamO.correctBoundaryConditions();
        lamN.correctBoundaryConditions();

        const scalarField& loI = lamO.primitiveField();
        const scalarField& liI = lamN.primitiveField();

        {
            scalarField& FgI = Fg.primitiveFieldRef();
            forAll(FgI, facei)
            {
                FgI[facei] *=
                    (FgI[facei] > 0)
                  ? Foam::min(loI[own[facei]], liI[nei[facei]])
                  : Foam::min(loI[nei[facei]], liI[own[facei]]);
            }
        }

        auto& FgBf = Fg.boundaryFieldRef();
        forAll(FgBf, patchi)
        {
            fvsPatchScalarField& pF = FgBf[patchi];
            const labelUList& fc = pF.patch().faceCells();

            if (pF.patch().coupled())
            {
                const scalarField loNei
                (
                    lamO.boundaryField()[patchi].patchNeighbourField()
                );
                const scalarField liNei
                (
                    lamN.boundaryField()[patchi].patchNeighbourField()
                );
                forAll(pF, bfacei)
                {
                    pF[bfacei] *=
                        (pF[bfacei] > 0)
                      ? Foam::min(loI[fc[bfacei]], liNei[bfacei])
                      : Foam::min(loNei[bfacei], liI[fc[bfacei]]);
                }
            }
            else
            {
                // physical patch: exports keep the positivity factor, INFLOW
                // is never throttled.. BCs outrank the limiter, throttling
                // BC inflow = volume without mass = the boundary air jet
                forAll(pF, bfacei)
                {
                    if (pF[bfacei] > 0)
                    {
                        pF[bfacei] *= loI[fc[bfacei]];
                    }
                }
            }
        }
    }

    // engagement counts for the log
    {
        label nOutC = 0, nInC = 0, nOutG = 0, nInG = 0;
        scalar mOutC = 1, mInC = 1, mOutG = 1, mInG = 1;
        const scalarField& oc = lamOutCond.primitiveField();
        const scalarField& ic = lamInCond.primitiveField();
        const scalarField& og = lamOutGas.primitiveField();
        const scalarField& ig = lamInGas.primitiveField();
        forAll(oc, celli)
        {
            if (oc[celli] < 1.0) { ++nOutC; mOutC = Foam::min(mOutC, oc[celli]); }
            if (ic[celli] < 1.0) { ++nInC;  mInC  = Foam::min(mInC,  ic[celli]); }
            if (og[celli] < 1.0) { ++nOutG; mOutG = Foam::min(mOutG, og[celli]); }
            if (ig[celli] < 1.0) { ++nInG;  mInG  = Foam::min(mInG,  ig[celli]); }
        }
        reduce(nOutC, sumOp<label>()); reduce(mOutC, minOp<scalar>());
        reduce(nInC,  sumOp<label>()); reduce(mInC,  minOp<scalar>());
        reduce(nOutG, sumOp<label>()); reduce(mOutG, minOp<scalar>());
        reduce(nInG,  sumOp<label>()); reduce(mInG,  minOp<scalar>());
        if (nOutC || nInC || nOutG || nInG)
        {
            Info<< "group flux limiter: condensed out " << nOutC
                << " (min " << mOutC << ") / in " << nInC
                << " (min " << mInC << "); gas out " << nOutG
                << " (min " << mOutG << ") / in " << nInG
                << " (min " << mInG << ")" << endl;
        }
    }

    // --- mass fractions  

    const dimensionedScalar cTinyY("cTinyY", dimDensity, 1e-16);

    // --- per species share of the limited group flux + fickian, then the
        //     conservative update. nothing else is alowed to touch c_i
    label nSpeciesLam = 0;
    {
        label i = 0;
        for (phaseModel& ph : phases_)
        {
            const surfaceScalarField& FgrpL = isGas[i] ? FGas : FCond;
            const volScalarField& cTot     = isGas[i] ? cGasTot : cCondTot;

            const volScalarField Yi
            (
                "Y." + ph.name(),
                max(c_[i], dimensionedScalar("z", dimDensity, 0.0))
               /max(cTot, cTinyY)
            );

            // the share: Y upwinded by the group flux sign, times it.
            // summed over the group this returns FgrpL exactly
            surfaceScalarField Fi
            (
                "F" + ph.name(),
                upwind<scalar>(mesh_, FgrpL).flux(Yi)
            );
            donorOutflow(Fi, FgrpL, Yi);    // outflow patch faces split by the cell's ledger, not the BC

            //  fickian diffusion within each group, listed pairs only..
            
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

            // per-species positivity: every export route (share, Fickian, patch) bounded by what the cell has
            {
                const scalarField& V   = mesh_.V().field();
                const labelUList& own  = mesh_.owner();
                const labelUList& nei  = mesh_.neighbour();
                const scalarField& cI  = c_[i].primitiveField();
                const scalarField& SuI = SuMass[i].primitiveField();

                scalarField outflow(mesh_.nCells(), 0.0);
                {
                    const scalarField& FiI = Fi.primitiveField();
                    forAll(FiI, facei)
                    {
                        if (FiI[facei] > 0) outflow[own[facei]] += FiI[facei];
                        else                outflow[nei[facei]] -= FiI[facei];
                    }
                    forAll(Fi.boundaryField(), patchi)
                    {
                        const fvsPatchScalarField& pF = Fi.boundaryField()[patchi];
                        const labelUList& fc = pF.patch().faceCells();
                        forAll(pF, bfacei)
                        {
                            if (pF[bfacei] > 0) outflow[fc[bfacei]] += pF[bfacei];
                        }
                    }
                }

                volScalarField lamSp
                (
                    IOobject("lam." + ph.name(), runTime.timeName(), mesh_),
                    mesh_, dimensionedScalar("one", dimless, 1.0)
                );
                scalarField& lamI = lamSp.primitiveFieldRef();
                forAll(lamI, celli)
                {
                    if (outflow[celli] > VSMALL)
                    {
                        const scalar budget =
                            Foam::max(cI[celli], scalar(0))*V[celli]*rDeltaT
                          + Foam::min(SuI[celli], scalar(0))*V[celli];
                        const scalar lam =
                            Foam::min(scalar(1), Foam::max(budget, scalar(0))/outflow[celli]);
                        if (lam < 1) { lamI[celli] = lam; ++nSpeciesLam; }
                    }
                }
                lamSp.correctBoundaryConditions();

                scalarField& FiI = Fi.primitiveFieldRef();
                forAll(FiI, facei)
                {
                    FiI[facei] *= (FiI[facei] > 0) ? lamI[own[facei]] : lamI[nei[facei]];
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
                            lamSp.boundaryField()[patchi].patchNeighbourField()
                        );
                        forAll(pF, bfacei)
                        {
                            pF[bfacei] *= (pF[bfacei] > 0) ? lamI[fc[bfacei]] : lamNei[bfacei];
                        }
                    }
                    else
                    {
                        forAll(pF, bfacei)
                        {
                            if (pF[bfacei] > 0) pF[bfacei] *= lamI[fc[bfacei]];
                        }
                    }
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

            // donor Cv_i and Cv_i*T^n ride with F_i, so TEqn moves the energy the ledger moved
            const volScalarField CvPh(ph.thermo().Cv());
            const volScalarField CvTPh(CvPh*T_);
            rhoCvPhi_ += Fi*donorValue(Fi, CvPh);
            rhoCvTPhi_ += Fi*donorValue(Fi, CvTPh);
            ++i;
        }
    }
    reduce(nSpeciesLam, sumOp<label>());
    if (nSpeciesLam)
    {
        Info<< "species positivity limiter: " << nSpeciesLam
            << " cell-species export(s) throttled" << endl;
    }


    //  now recover alphas from the mass-balance. c_i can carry a tiny negative
    //  from roundoff hence the max(..,0) everywhere

 

    // purge
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


    // hierarchical recovery of the alphas (a global rescale failed): alpha = c/rho sums to 1 only if the
    // transported masses match the EOS densities at the current p and T, which they never quite do
    // (species move with last step's fluxes while p moves by bars), so a defect builds every step

    const dimensionedScalar rhoMinRec("rhoMinRec", dimDensity, 1e-3);

    // closureGasFracMin: the gas mass balance must explain at least this share of the volume it is asked to
    // fill before it is trusted (VSMALL made t/tiny huge); closureCondNearFull: above this condensed fill the
    // leftover is a splitting defect, not a void, and the condensed phase absorbs it
    const scalar gasFracMin =
        phasedictionary.getOrDefault<scalar>("closureGasFracMin", 0.1);
    const scalar condNearFull =
        phasedictionary.getOrDefault<scalar>("closureCondNearFull", 0.95);

    // below this much condense phase a cell counts as a void 
    const scalar condVoidFloor = 1e-3;

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
                         && (vc >= condNearFull || vc >= Foam::max(vg, condVoidFloor))
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
                        else if (vc >= condNearFull || vc >= Foam::max(vg, condVoidFloor))
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
                // only count absorbs that mean something
                if
                (
                    Vg > 1e-6
                 && vGasI[celli] < gasFracMin*Vg
                 && (vc >= condNearFull || vc >= Foam::max(vGasI[celli], condVoidFloor))
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

    // refresh the mass-balance patch values, c_b = alpha_b*rho_b on physical
    // patches 


    {
        label i = 0;
        for (phaseModel& ph : phases_)
        {
            const volScalarField::Boundary& aBf = ph.boundaryField();
            const volScalarField& rhoPh = ph.thermo().rho();
            const volScalarField::Boundary& rBf = rhoPh.boundaryField();

            volScalarField::Boundary& cBf = c_[i].boundaryFieldRef();

            forAll(cBf, patchi)
            {
                if (cBf[patchi].coupled())
                {
                    // halo, gets the proper swap below
                    continue;
                }
                cBf[patchi] == aBf[patchi]*rBf[patchi];
            }

            // true halo exchange on the processor/cyclic patches
            c_[i].correctBoundaryConditions();
            ++i;
        }
    }

    // the volume closure is evaluated after TEqn (volumeClosure), on the state the
    // pressure equation actually starts from; PCR returned here is phase change only
    {
        static bool warnedOldClosure = false;
        if (!warnedOldClosure)
        {
            warnedOldClosure = true;
            for
            (
                const word key
              : {"closureRelax", "closureCoFraction", "closureVolLimit", "closureMaskPatches"}
            )
            {
                if (phasedictionary.found(key))
                {
                    WarningInFunction
                        << key << " is ignored: the volume closure is evaluated after"
                        << " TEqn, unrelaxed, uncapped and unmasked" << endl;
                }
            }
        }
    }

    // ledger trace at (p^n, T^n): worst under-/over-full cell and everything moving ledger mass there
    {
        const scalar traceBand =
            phasedictionary.getOrDefault<scalar>("closureTraceBand", 0.1);

        const scalarField& sI = sumImplied.primitiveField();
        label iMin = -1, iMax = -1;
        scalar sMin = GREAT, sMax = -GREAT;
        forAll(sI, celli)
        {
            if (sI[celli] < sMin) { sMin = sI[celli]; iMin = celli; }
            if (sI[celli] > sMax) { sMax = sI[celli]; iMax = celli; }
        }
        scalar gsMin = sMin, gsMax = sMax;
        reduce(gsMin, minOp<scalar>());
        reduce(gsMax, maxOp<scalar>());

        const bool traceMin = (gsMin < 1.0 - traceBand);
        const bool traceMax = (gsMax > 1.0 + traceBand);

        if (traceMin || traceMax)
        {
            // per unit cell volume, per step: positive = leaves the cell
            const scalarField divPhiDt(dt*fvc::div(phi_)().primitiveField());
            const scalarField divFGasDt(dt*fvc::div(FGas)().primitiveField());
            const scalarField divFCondDt(dt*fvc::div(FCond)().primitiveField());

            auto trace = [&](const word& tag, const label celli)
            {
                scalar betaRec = 0, suGasDt = 0, suCondDt = 0;
                label i = 0;
                for (const phaseModel& ph : phases_)
                {
                    if (isGas[i]) suGasDt += dt*SuMass[i][celli];
                    else { betaRec += ph[celli]; suCondDt += dt*SuMass[i][celli]; }
                    ++i;
                }

                Pout<< "ledger trace (p^n, T^n) [" << tag << "] cell " << celli
                    << " at " << mesh_.C()[celli]
                    << " s=" << sI[celli]
                    << " vCond=" << vCond[celli] << " vGas=" << vGas[celli]
                    << " beta0=" << condensate[celli] << " betaRec=" << betaRec
                    << " aI=" << aI[celli]
                    << " T=" << T_[celli] << " p=" << p_[celli] << nl
                    << "    per step, per cell volume, +ve = out:"
                    << " dt*div(phi)=" << divPhiDt[celli]
                    << " dt*div(FGas)=" << divFGasDt[celli]
                    << " dt*div(FCond)=" << divFCondDt[celli]
                    << " dt*Su(gas)=" << suGasDt
                    << " dt*Su(cond)=" << suCondDt << endl;

                i = 0;
                for (const phaseModel& ph : phases_)
                {
                    const scalar rhoEOS = ph.thermo().rho()()[celli];
                    const scalar cI = c_[i][celli];
                    Pout<< "    " << ph.name()
                        << ": c=" << cI
                        << " alphaRec=" << ph[celli]
                        << " rhoEOS=" << rhoEOS
                        << " c/alpha=" << cI/Foam::max(ph[celli], scalar(1e-12))
                        << " v=c/rho=" << cI/Foam::max(rhoEOS, rhoMinRec.value())
                        << " dt*Su=" << dt*SuMass[i][celli] << endl;
                    ++i;
                }
            };

            if (traceMin && iMin >= 0 && sMin == gsMin) trace("under-full", iMin);
            if (traceMax && iMax >= 0 && sMax == gsMax) trace("over-full", iMax);
        }
    }


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


Foam::tmp<Foam::volScalarField>
Foam::multiphaseMixtureThermo::volumeClosure(volScalarField& chi) const
{
    const Time& runTime = mesh_.time();
    const dimensionedScalar rhoMinRec("rhoMinRec", dimDensity, 1e-3);
    const dimensionedScalar cZero("cZero", dimDensity, Zero);

    // s: volume the ledger needs at the current (p, T), per unit cell volume
    volScalarField s
    (
        IOobject("closureS", runTime.timeName(), mesh_),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    // rhoCv: what the p*dV work heats; pdsdT = p*ds/dT|p: the volume that heat buys back
    volScalarField rhoCvL
    (
        IOobject("closureRhoCv", runTime.timeName(), mesh_),
        mesh_,
        dimensionedScalar(dimEnergy/dimVolume/dimTemperature, Zero)
    );
    volScalarField pdsdT
    (
        IOobject("closurePdsdT", runTime.timeName(), mesh_),
        mesh_,
        dimensionedScalar(dimEnergy/dimVolume/dimTemperature, Zero)
    );

    label i = 0;
    for (const phaseModel& ph : phases_)
    {
        const rhoThermo& thermo = ph.thermo();
        const volScalarField cPos(max(c_[i], cZero));
        const volScalarField vi(cPos/max(thermo.rho(), rhoMinRec));

        s += vi;

        const volScalarField CvPh(thermo.Cv());
        rhoCvL += cPos*CvPh;

        // isobaric expansivity from Cp - Cv = T*beta^2/(rho*kappaT), kappaT = psi/rho:
        // 1/T for perfectGas/perfectFluid(rho0 0), 0 for adiabaticPerfectFluid
        const volScalarField CpMCv(thermo.Cp() - CvPh);
        const volScalarField betaPh
        (
            sqrt
            (
                max(CpMCv, dimensionedScalar(CpMCv.dimensions(), Zero))
               *thermo.psi()/T_
            )
        );
        pdsdT += p_*vi*betaPh;

        ++i;
    }

    // small capacity so an emptied cell refills with chi -> 1 instead of 0/0
    const dimensionedScalar rhoCvSmall(rhoCvL.dimensions(), 1e-3);
    chi = (rhoCvL + rhoCvSmall)/(rhoCvL + pdsdT + rhoCvSmall);
    chi.correctBoundaryConditions();

    const dimensionedScalar rDeltaT(1.0/runTime.deltaT());

    tmp<volScalarField> tSrc
    (
        new volScalarField
        (
            IOobject("volumeClosureSrc", runTime.timeName(), mesh_),
            chi*(s - scalar(1))*rDeltaT
        )
    );

    const scalarField& sI = s.primitiveField();
    const scalarField& chiI = chi.primitiveField();
    Info<< "volume closure (post-TEqn): s min/max = " << gMin(sI)
        << " / " << gMax(sI)
        << ", chi min/max = " << gMin(chiI) << " / " << gMax(chiI)
        << ", max|chi*(s - 1)| = "
        << gMax(mag(chiI*(sI - 1.0))()) << endl;

    return tSrc;
}


// ************************************************************************* //