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

#include "fvmSup.H"


#include "fvcDdt.H"
#include "fvcFlux.H"
#include "fvcMeshPhi.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
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


void Foam::multiphaseMixtureThermo::buildPhaseIndexMap()
{
    label idx = 0;
    for (const phaseModel& phase : phases_)
    {
        phaseIndexMap_.insert(phase.name(), idx);
        ++idx;
    }
}


Foam::label Foam::multiphaseMixtureThermo::phaseIndex
(
    const word& phaseName
) const
{
    const auto iter = phaseIndexMap_.cfind(phaseName);

    if (!iter.good())
    {
        FatalErrorInFunction
            << "Phase '" << phaseName << "' not found in phase list."
            << nl << "Available phases: " << phaseIndexMap_.sortedToc()
            << exit(FatalError);
    }

    return iter();
}


void Foam::multiphaseMixtureThermo::readPhaseTransitions()
{
    if (!found("phaseTransitions"))
    {
        Info<< "No phaseTransitions sub-dictionary found. "
            << "No state transitions will be modelled." << nl << endl;
        return;
    }

    const dictionary& transDict = subDict("phaseTransitions");
    wordList transNames(transDict.toc());

    phaseTransitions_.resize(transNames.size());

    label ti = 0;
    for (const word& tName : transNames)
    {
        const dictionary& td = transDict.subDict(tName);

        phaseTransition& pt = phaseTransitions_[ti];
        pt.name = tName;

        const word typeStr(td.get<word>("type"));

        if (typeStr == "boiling")
        {
            pt.type = transitionType::BOILING;
            pt.phase1Name = td.get<word>("condensedPhase");
            pt.phase2Name = td.get<word>("gaseousPhase");
            pt.boilingTemperature = td.get<scalar>("boilingTemperature");
            pt.latentHeat = td.get<scalar>("latentHeat");

            // Validate that both phases exist
            phaseIndex(pt.phase1Name);
            phaseIndex(pt.phase2Name);

            Info<< "Registered phase transition '" << tName << "': "
                << "boiling  " << pt.phase1Name << " <-> " << pt.phase2Name
                << "  T_boil = " << pt.boilingTemperature << " K"
                << "  L = " << pt.latentHeat << " J/kg"
                << endl;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown transition type '" << typeStr
                << "' in phaseTransitions/" << tName << nl
                << "Valid types: boiling"
                << exit(FatalError);
        }

        ++ti;
    }

    Info<< "Registered " << phaseTransitions_.size()
        << " phase transition(s)." << nl << endl;
}


void Foam::multiphaseMixtureThermo::initAlphaRho()
{
    alphaRho_.resize(phases_.size());

    label idx = 0;
    for (const phaseModel& phase : phases_)
    {
        // Build BC types: use fixedValue for physical boundaries,
        // keep special geometric types (empty, wedge, processor, etc).
        const volScalarField::Boundary& alphaBf = phase.boundaryField();
        wordList alphaRhoBCTypes(alphaBf.size());

        forAll(alphaBf, patchi)
        {
            const word& bcType = alphaBf[patchi].type();

            if
            (
                bcType == "empty"
             || bcType == "wedge"
             || bcType == "symmetry"
             || bcType == "symmetryPlane"
             || bcType == "processor"
             || bcType == "processorCyclic"
             || bcType == "cyclic"
             || bcType == "cyclicAMI"
            )
            {
                alphaRhoBCTypes[patchi] = bcType;
            }
            else
            {
                // Walls, inlets, outlets, alphaContactAngle, etc.
                // → fixedValue so the implicit solver respects
                //   the values we set from alpha * rho
                alphaRhoBCTypes[patchi] = "fixedValue";
            }
        }

        alphaRho_.set
        (
            idx,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("alphaRho", phase.name()),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimDensity, 0.0),
                alphaRhoBCTypes
            )
        );

        // Set internal field
        alphaRho_[idx].primitiveFieldRef() =
            phase.primitiveField()
          * phase.thermo().rho()().primitiveField();

        // Set boundary face values from alpha * rho
        {
            const tmp<volScalarField> trhoK(phase.thermo().rho());
            const volScalarField& rhoK = trhoK();

            forAll(alphaRho_[idx].boundaryFieldRef(), patchi)
            {
                if (alphaRhoBCTypes[patchi] == "fixedValue")
                {
                    alphaRho_[idx].boundaryFieldRef()[patchi] ==
                        phase.boundaryField()[patchi]
                      * rhoK.boundaryField()[patchi];
                }
            }
        }

        alphaRho_[idx].correctBoundaryConditions();
        alphaRho_[idx].oldTime();

        // Handle restarting the solver and maintaining the previous values
        // (only needs to run if solver has restarted from latestTime.)
        if (bFreshStart && mesh_.time().timeIndex() > 1)
        {
            volScalarField alphaRhoOldTime_idx_
            (
                IOobject
                (
                    "alphaRho." + phase.name() + ".oldTime",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                alphaRho_[idx]
            );

            if (alphaRhoOldTime_idx_.headerOk())
            {
                Info<< "Mixture: Restoring oldTime values "
                    << "for alphaRho." << phase.name()
                    << endl;

                alphaRho_[idx].oldTime() = alphaRhoOldTime_idx_;
            }
            else
            {
                Info<< "Mixture: WARNING could not find field "
                    << "for `alphaRho." << phase.name() << ".oldTime`. "
                    << "This will cause issues in fvc::ddt terms!"
                    << "PLEASE check your latestTime directory (folder)."
                    << endl;
            }

            // bFreshStart is disabled after this loop!
        }

        Info<< "Initialised alphaRho." << phase.name()
            << "  mass = "
            << gSum(alphaRho_[idx].primitiveField() * mesh_.V().field())
            << " kg" << endl;

        ++idx;
    }

    if (bFreshStart)
    {
        // close it off after the first read through
        bFreshStart = false;
    }
}

void Foam::multiphaseMixtureThermo::writeOldTimeValues(void)
{
    // Update any old time storage values related to `.oldTime()` calls
    // Presently, this is just for `alphaRho_`.
    if (mesh_.time().writeTime())
    {
        label idx = 0;
        for (const phaseModel& phase : phases_)
        {
            // alphaRho_[idx].correctBoundaryConditions();

            volScalarField alphaRhoOldTime_idx_
            (
                IOobject
                (
                    "alphaRho." + phase.name() + ".oldTime",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                alphaRho_[idx]
            );
            // alphaRhoOldTime_idx_.correctBoundaryConditions();
            alphaRhoOldTime_idx_ = alphaRho_[idx].oldTime();
            alphaRhoOldTime_idx_.write();

            idx++;
        }
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
            IOobject::READ_IF_PRESENT, //NO_READ,
            IOobject::AUTO_WRITE //NO_WRITE
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

    buildPhaseIndexMap();
    readPhaseTransitions();
    initAlphaRho();

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


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::rhoAlphaRho() const
{
    tmp<volScalarField> trho
    (
        new volScalarField
        (
            IOobject
            (
                "rhoAlphaRho",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("zero", dimDensity, 0.0)
        )
    );

    volScalarField& rhoField = trho.ref();

    forAll(alphaRho_, idx)
    {
        rhoField += alphaRho_[idx];
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
    volScalarField *massdotterm,
    bool finalIter
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
            PCR = solveAlphas(massdotterm, finalIter);
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi_;
        }

        rhoPhi_ = rhoPhiSum;
    }
    else
    {
        PCR = solveAlphas(massdotterm, finalIter);
    }

    // write all `alphaRho.phase.oldTime()` to file for restart capability.
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

    return gradAlphaf/(mag(gradAlphaf) + deltaN_);
}


Foam::tmp<Foam::surfaceScalarField> Foam::multiphaseMixtureThermo::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}


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

            if (uTheta > SMALL)
            {
                const scalar thetaA = degToRad(tp().thetaA(matched));
                const scalar thetaR = degToRad(tp().thetaR(matched));

                vectorField Uwall
                (
                    U_.boundaryField()[patchi].patchInternalField()
                  - U_.boundaryField()[patchi]
                );
                Uwall -= (AfHatPatch & Uwall)*AfHatPatch;

                vectorField nWall
                (
                    nHatPatch - (AfHatPatch & nHatPatch)*AfHatPatch
                );

                nWall /= (mag(nWall) + SMALL);

                scalarField uwall(nWall & Uwall);

                theta += (thetaA - thetaR)*tanh(uwall/uTheta);
            }

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


// void Foam::multiphaseMixtureThermo::updateRhoPhi()
// {
//     rhoPhi_ = dimensionedScalar(dimensionSet(1, 0, -1, 0, 0), Zero);

//     label idx = 0;
//     for (const phaseModel& phase : phases_)
//     {
//         rhoPhi_ += fvc::interpolate(alphaRho_[idx]) * phi_;
//         ++idx;
//     }
// }

// void Foam::multiphaseMixtureThermo::syncAlphaRho()
// {
//     label idx = 0;
//     for (const phaseModel& phase : phases_)
//     {
//         const volScalarField& rhoK = phase.thermo().rho();

//         alphaRho_[idx].primitiveFieldRef() =
//             phase.primitiveField() * rhoK.primitiveField();

//         forAll(alphaRho_[idx].boundaryFieldRef(), patchi)
//         {
//             if (alphaRho_[idx].boundaryField()[patchi].type() == "fixedValue")
//             {
//                 alphaRho_[idx].boundaryFieldRef()[patchi] ==
//                     phase.boundaryField()[patchi]
//                   * rhoK.boundaryField()[patchi];
//             }
//         }

//         alphaRho_[idx].correctBoundaryConditions();

//         ++idx;
//     }
// }

Foam::tmp<Foam::volScalarField>
Foam::multiphaseMixtureThermo::continuityError() const
{
    tmp<volScalarField> tErr
    (
        new volScalarField
        (
            IOobject
            (
                "continuityError",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(dimDensity/dimTime, Zero)
        )
    );

    volScalarField& err = tErr.ref();

    forAll(alphaRho_, idx)
    {
        err += fvc::ddt(alphaRho_[idx]);
    }

    err += fvc::div(rhoPhi_);

    return tErr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
//  solveAlphas — implicit conservative partial-density transport
//
//  Transports alphaRho_k = alpha_k * rho_k implicitly:
//
//      fvm::ddt(alphaRho_k) + fvm::div(phi, alphaRho_k) = Mdot_k ...... with diffusion and compression terms now re-added
//
//
//  The implicit solver handles compressibility via the time
//  derivative.
//
//  rhoPhi is extracted from the solved matrix via .flux(), which
//  gives the face fluxes exactly consistent with the field update.
//  This guarantees:
//      d(rho)/dt + div(rhoPhi) = 0 ! ! ! ! ! ! !
//  so fvm::ddt(rho, T) + fvm::div(rhoPhi, T) has zero spurious
//  source contributions.  Temperature/energy solve stability is guaranteed by construction.
//
//  Per-pair mass is conserved to solver tolerance/ machine precision.
//  No (divU-PCR)*alpha is required in this formulation - although looking back I think I fucked
// that up anyway as the compressability was being applied inconsistently to all componnets in cells
//
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::multiphaseMixtureThermo::solveAlphas
(
    volScalarField *massdotterm,
    bool finalIter
)
{
    tmp<volScalarField> tPCR
    (
        new volScalarField
        (
            IOobject("PhaseChangeRate", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedScalar("PhaseChangeRate", dimensionSet(0,0,-1,0,0), 0.0)
        )
    );
    volScalarField& PCR = tPCR.ref();

    *massdotterm *= 0.0;

    const label nPhases = phases_.size();


    //need condensate for Hertz-Knudsen as rate depends on nucleation sites
    volScalarField condensate
    (
        IOobject("condensate", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar("condensate", dimless, 1.0)
    );

    for (const phaseModel& alpha : phases_)
    {
        if (alpha.isGaseous())
        {
            condensate -= alpha;
        }
    }




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

    const dimensionedScalar P0("P0", dimensionSet(1,-1,-2,0,0), phasedictionary);
    const dimensionedScalar int_thickness("interface_thickness", dimensionSet(0,1,0,0,0), phasedictionary);
    const dimensionedScalar gasconstant("gasconstant", dimensionSet(1,2,-2,-1,-1), scalar(8.314));
    const dimensionedScalar maxrate("maxrate", dimensionSet(0,0,-1,0,0), 0.5/mesh_.time().deltaT().value());

    PtrList<volScalarField> Mdot(nPhases);
    {
        label idx = 0;
        for (const phaseModel& phase : phases_)
        {
            Mdot.set
            (
                idx,
                new volScalarField
                (
                    IOobject("Mdot_" + phase.name(), mesh_.time().timeName(), mesh_),
                    mesh_,
                    dimensionedScalar("zero", dimDensity/dimTime, 0.0)
                )
            );
            ++idx;
        }
    }

    for (const phaseTransition& trans : phaseTransitions_)
    {
        switch (trans.type)
        {
            case transitionType::BOILING:
            {
                const label idxLiq = phaseIndex(trans.phase1Name);
                const label idxVap = phaseIndex(trans.phase2Name);

                const phaseModel& phaseLiq = phases_[trans.phase1Name];
                const phaseModel& phaseVap = phases_[trans.phase2Name];

                const volScalarField& alphaLiq = phaseLiq;
                const volScalarField& alphaVap = phaseVap;

                const volScalarField rhoLiq(phaseLiq.thermo().rho());
                const volScalarField rhoVap(phaseVap.thermo().rho());

                const dimensionedScalar pair_boil_T("pair_boil_T", dimensionSet(0,0,0,1,0), trans.boilingTemperature);
                const dimensionedScalar pair_LHG("pair_LHG", dimensionSet(0,2,-2,0,0), trans.latentHeat);

                volScalarField Psat
                (
                    IOobject("Psat", mesh_.time().timeName(), mesh_),
                    mesh_,
                    dimensionedScalar("Psat", dimensionSet(1,-1,-2,0,0), 0.0)
                );
                Psat = P0 * Foam::exp(((phaseLiq.thermo().W()/1000.0)*pair_LHG/(pair_boil_T*gasconstant))*(1.0 - (pair_boil_T/T_)));///1000 as standard openfoam units are for some reason g/mol and we want kg/mol

                volScalarField evapcoefffield
                (
                    IOobject("evapcoefffield", mesh_.time().timeName(), mesh_),
                    mesh_,
                    dimensionedScalar("zero", dimless, 0.0),
                    zeroGradientFvPatchScalarField::typeName
                );
                forAll(mesh_.C(), celli)
                {
                    evapcoefffield[celli] = (p_[celli] > Psat[celli]) ? 0.0 : 1.0;
                }

                const volScalarField evaprate
                (
                    Foam::sqrt((phaseLiq.thermo().W()/1000.0)/(2.0*M_PI*gasconstant*T_))
                  * (1.0/(int_thickness*rhoLiq)) * (Psat - p_)
                );

                const volScalarField condrate
                (
                    -max(condensate, 1e-6)
                  * Foam::sqrt((phaseVap.thermo().W()/1000.0)/(2.0*M_PI*gasconstant*T_))
                  * (1.0/(int_thickness*rhoVap)) * (Psat - p_)
                );

                const volScalarField evapContrib(min(evaprate, maxrate) * evapcoefffield * alphaLiq);
                const volScalarField condContrib(min(condrate, maxrate) * (1.0 - evapcoefffield) * alphaVap);
                const volScalarField densityRatio(rhoVap / rhoLiq);

                // Single mass rate — applied symmetrically! I think this should be quicker than the double symmetrical loop thing I had before
                const volScalarField netMassRate(rhoVap * (evapContrib - condContrib));

                Mdot[idxLiq] -= netMassRate;
                Mdot[idxVap] += netMassRate;

                // Volumetric dilation passed for pressure equation....
                PCR += netMassRate/rhoVap - netMassRate/rhoLiq;

                // Energy source passed to the energy equation eventually.....
                const volScalarField CvLiq(phaseLiq.thermo().Cv());
                const volScalarField CvVap(phaseVap.thermo().Cv());

                *massdotterm += pair_LHG *
                (
                    (rhoLiq/CvLiq)*condContrib*densityRatio
                  - (rhoVap/CvVap)*evapContrib*densityRatio
                );
                *massdotterm += pair_LHG *
                (
                    (rhoVap/CvVap)*evapContrib
                  - (rhoLiq/CvLiq)*condContrib
                );

                Info<< "Boiling transition: ("
                    << trans.phase1Name << " <-> " << trans.phase2Name
                    << ")  T_boil = " << pair_boil_T.value() << " K" << endl;
                break;
            }
        }
    }

    // Interface compression & diffusion
    PtrList<surfaceScalarField> phiR(nPhases);

    const bool bEnableInterfaceCompression = true;
    const scalar cAlphaMax = 2.0;

    if (bEnableInterfaceCompression)
    {
        label idx = 0;
        for (const phaseModel& alpha : phases_)
        {
            phiR.set
            (
                idx,
                new surfaceScalarField
                (
                    IOobject
                    (
                        "phiR." + alpha.name(),
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar(dimVol/dimTime, Zero)
                )
            );

            const volScalarField rhoK(alpha.thermo().rho());

            for (const phaseModel& alpha2 : phases_)
            {
                if (&alpha2 == &alpha) continue;

                // --- Interface compression ---
                auto cAlpha = cAlphas_.cfind(interfacePair(alpha, alpha2));
                if (cAlpha.good())
                {
                    // Info<< "Debug: cAlpha = " << scalar(cAlpha()) << endl;
                    surfaceScalarField phic
                    (
                        mag(phi_) / mesh_.magSf()
                    );

                    phiR[idx] += min
                    (
                        scalar(cAlpha()) * phic,
                        cAlphaMax * max(phic)
                    ) * nHatf(alpha, alpha2)
                    ;
                }

                // --- Interface diffusion ---
                auto dAlpha = dAlphas_.cfind(interfacePair(alpha, alpha2));
                if (dAlpha.good())
                {
                    dimensionedScalar valdiff("valdiff", dimdiff_, dAlpha());

                    phiR[idx] -=
                    (
                          valdiff
                        * mesh_.magSf()
                        * (
                              fvc::interpolate(alpha2) * fvc::snGrad(alpha)
                            - fvc::interpolate(alpha)  * fvc::snGrad(alpha2)
                        )
                    );
                }
            }

            ++idx;
        }
    }


    // ================================================================
    //   Implicit transport of alphaRho with flux extraction for the msin solver transport equations
    //
    //  fvm::ddt(alphaRho) + fvm::div(phi, alphaRho) = Mdot
    //
    //  The matrix flux (alphaRhoEqn.flux()) is the face flux that
    //  is exactly consistent with the implicit field update.
    //
    //  rhoPhi = SUM(alphaRhoEqn.flux())
    //
    //  This guarantees d(rho)/dt + div(rhoPhi) = 0 to solver
    //  tolerance, so TEqn has no spurious density source.
    //
    //  Note: with interface compression enabled, a mass source can
    //  emerge (see eq. below), this can break conservation slightly.
    //
    //      d(rho)/dt + div(rhoPhi) = -div(rhoPhiR)
    //   => d(rho_k alpha_k)/dt + div(rho_k alpha_k U) =
    //                                         Mdot[k] - div(rho_k alpha_k U_R)
    //
    //  where the latter is the phase-wise continuity, which includes Mdot from
    //  the state
    // ================================================================

    rhoPhi_ = dimensionedScalar(dimensionSet(1, 0, -1, 0, 0), Zero);

    {
        label idx = 0;
        for (phaseModel& phase : phases_)
        {
            // Sync boundary values of alphaRho from alpha * rho.
            // Only fixedValue patches — empty/processor/wedge handle
            // themselves through correctBoundaryConditions().
            forAll(alphaRho_[idx].boundaryFieldRef(), patchi)
            {
                if (alphaRho_[idx].boundaryField()[patchi].type() == "fixedValue")
                {
                    alphaRho_[idx].boundaryFieldRef()[patchi] ==
                        phase.boundaryField()[patchi]
                      * phase.thermo().rho().boundaryField()[patchi];
                }
            }

            // const volScalarField rhoK(phase.thermo().rho());

fvScalarMatrix alphaRhoEqn
(
    fvm::ddt(alphaRho_[idx])
  + fvm::div(phi_, alphaRho_[idx], "div(phi,alphaRho)")
  + fvm::div(phiR[idx], alphaRho_[idx], "div(phiR,alphaRho)")
 ==
    Mdot[idx]
);

            alphaRhoEqn.solve();

            // Extract the EXACT face flux from the solved matrix
            rhoPhi_ += alphaRhoEqn.flux();

            Info<< "Solved alphaRho." << phase.name()
                << "  mass = "
                << gSum(alphaRho_[idx].primitiveField() * mesh_.V().field())
                << " kg" << endl;

            ++idx;
        }
    }


    // ================================================================
    //  Recover alpha = alphaRho / rho.  NO renormalisation.
    //
    //  Sigma(alpha) != 1 is the PHYSICAL volume deficit/surplus.
    //  When vapour condenses, Sigma(alpha) < 1 because the liquid
    //  occupies less volume.  This deficit is passed to the pressure
    //  equation via vDot (PCR), which adjusts div(U) to drive air
    //  inflow from the boundaries.  Over PIMPLE iterations, phi
    //  converges and Sigma(alpha) approaches 1.

    // ================================================================

    volScalarField sumAlpha
    (
        IOobject("sumAlpha", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    {
        label idx = 0;
        for (phaseModel& phase : phases_)
        {
            // const volScalarField& rhoK = phase.thermo().rho();
            const volScalarField& rhoK(phase.thermo().rho());

            phase.primitiveFieldRef() =
                max
                (
                    alphaRho_[idx].primitiveField() / rhoK.primitiveField(),
                    scalar(0)
                );

            phase.correctBoundaryConditions();

            Info<< phase.name() << " volume fraction, min, max = "
                << phase.weightedAverage(mesh_.V()).value()
                << ' ' << min(phase).value()
                << ' ' << max(phase).value()
                << endl;

            sumAlpha += phase;
            ++idx;
        }
    }

    Info<< "Phase-sum volume fraction, min, max = "
        << sumAlpha.weightedAverage(mesh_.V()).value()
        << ' ' << min(sumAlpha).value()
        << ' ' << max(sumAlpha).value()
        << endl;



        // ================================================================
// Operator-split compression/diffusion on alpha
//
//  Compression is a volume concept (sharpens alpha profile).
//   cannot be embedded in independent per-phase alphaRho
//  equations because  no cross-phase coupling to
//  enforce sum(alpha) = 1 (MULES uses limitSum for this in the old formulation).
//
//  Applied as bounded explicit correction to alpha, then
//  alphaRho is rebuilt.  The mass change from rebuilding is
//  physical: it corresponds to interface sharpening between
//  materials of different density, and vanishes with mesh
//  refinement.  rhoPhi is updated for TEqn consistency.
// ================================================================
//////////////////////////// THIS SECTION IS FOR DIFFUSION AND COMPRESSION
// if (finalIter)
// {

// {
//     const scalar dt = mesh_.time().deltaTValue();

//     // --- Build per-phase compression + diffusion alpha fluxes ---

//     PtrList<surfaceScalarField> alphaPhiR(nPhases);

//     {
//         label idx = 0;
//         for (const phaseModel& alpha : phases_)
//         {
//             // Compression velocity for this phase
//             surfaceScalarField phiRAlpha
//             (
//                 IOobject
//                 (
//                     "phiRAlpha." + alpha.name(),
//                     mesh_.time().timeName(),
//                     mesh_
//                 ),
//                 mesh_,
//                 dimensionedScalar(dimVol/dimTime, Zero)
//             );

//             for (const phaseModel& alpha2 : phases_)
//             {
//                 if (&alpha2 == &alpha) continue;

//                 auto cAlpha = cAlphas_.cfind(interfacePair(alpha, alpha2));
//                 if (cAlpha.good())
//                 {
//                     surfaceScalarField phic(mag(phi_) / mesh_.magSf());

//                     phiRAlpha += min
//                     (
//                         scalar(cAlpha()) * phic,
//                         max(phic)
//                     ) * nHatf(alpha, alpha2);
//                 }

//                 auto dAlpha = dAlphas_.cfind(interfacePair(alpha, alpha2));
//                 if (dAlpha.good())
//                 {
//                     dimensionedScalar valdiff("valdiff", dimdiff_, dAlpha());

//                     phiRAlpha -= valdiff
//                         * mesh_.magSf()
//                         * (
//                               fvc::interpolate(alpha2)
//                             * fvc::snGrad(alpha)
//                             - fvc::interpolate(alpha)
//                             * fvc::snGrad(alpha2)
//                           );
//                 }
//             }

//             // Face flux of alpha driven by compression/diffusion
//             alphaPhiR.set
//             (
//                 idx,
//                 new surfaceScalarField
//                 (
//                     "alphaPhiR." + alpha.name(),
//                     fvc::flux
//                     (
//                         phiRAlpha,
//                         static_cast<const volScalarField&>(alpha),
//                         "div(phiR,alpha)"
//                     )
//                   - fvc::interpolate
//                     (
//                         static_cast<const volScalarField&>(alpha)
//                     ) * phiRAlpha
//                 )
//             );

//             ++idx;
//         }
//     }

//     //  Apply bounded correction to alpha, rebuild alphaRho

//     {
//         label idx = 0;
//         for (phaseModel& phase : phases_)
//         {
//             volScalarField& alpha = phase;
//             const volScalarField rhoK(phase.thermo().rho());

//             volScalarField deltaAlpha
// (
//     mesh_.time().deltaT() * fvc::div(alphaPhiR[idx])
// );

//             // Bound so alpha stays in [0, 1]
//             deltaAlpha = max(deltaAlpha, -alpha);
//             deltaAlpha = min(deltaAlpha, scalar(1) - alpha);

//             alpha += deltaAlpha;
//             alpha.correctBoundaryConditions();

//             // Rebuild conserved field from corrected alpha
//             alphaRho_[idx].primitiveFieldRef() =
//                 alpha.primitiveField() * rhoK.primitiveField();

//             forAll(alphaRho_[idx].boundaryFieldRef(), patchi)
//             {
//                 if (alphaRho_[idx].boundaryField()[patchi].type()
//                     == "fixedValue")
//                 {
//                     alphaRho_[idx].boundaryFieldRef()[patchi] ==
//                         alpha.boundaryField()[patchi]
//                       * rhoK.boundaryField()[patchi];
//                 }
//             }
//             alphaRho_[idx].correctBoundaryConditions();

//             // Update rhoPhi for TEqn consistency
//             rhoPhi_ += fvc::interpolate(rhoK) * alphaPhiR[idx];

//             ++idx;
//         }
//     }
// }

// }
//////////////////////////// THIS SECTION IS FOR DIFFUSION AND COMPRESSION



    // ================================================================
    //  The deviation (Sigma(alpha) - 1) feeds into PCR:
    //    Sigma(alpha) < 1 (void) → PCR correction < 0
    //      → pEqn creates deeper pressure minimum at void
    //      → pressure gradient drives air inflow from boundary
    //    Sigma(alpha) > 1 (excess) → PCR correction > 0
    //      → pEqn creates pressure maximum → pushes excess out
    //
    //  Sign note: the incompressible code uses (1 - Σα) because
    //  there the correction enters the alpha equation's Su term
    //  as -PCR*alpha (opposite sign convention).  Here the
    //  correction goes directly to pEqn via vDot, so the sign
    //  is (Σα - 1).
    //
    //  This does NOT modify alphaRho — mass stays exactly conserved.
    // ================================================================

    const scalar alphaVolCorrCoeff = 1.0;

    PCR += (sumAlpha - 1.0) * alphaVolCorrCoeff
         / mesh_.time().deltaT();

    Info<< "Volume imbalance: mean(sumAlpha-1) = "
        << (sumAlpha - 1.0)().weightedAverage(mesh_.V()).value()
        << endl;


    // ================================================================
    // Mass diagnostics from alphaRho (the conserved field)
    // ================================================================

    Info << "Phase pair masses in domain (from alphaRho):" << endl;

    HashTable<bool> processedPhases;

    for (const phaseTransition& trans : phaseTransitions_)
    {
        const scalar mass1 =
            gSum(alphaRho_[phaseIndex(trans.phase1Name)].primitiveField()
               * mesh_.V().field());
        const scalar mass2 =
            gSum(alphaRho_[phaseIndex(trans.phase2Name)].primitiveField()
               * mesh_.V().field());

        Info << "    " << trans.phase1Name << " + " << trans.phase2Name
             << ": " << (mass1 + mass2) << " kg ("
             << trans.phase1Name << ": " << mass1 << " kg, "
             << trans.phase2Name << ": " << mass2 << " kg)"
             << endl;

        processedPhases.insert(trans.phase1Name, true);
        processedPhases.insert(trans.phase2Name, true);
    }

    for (const phaseModel& phase : phases_)
    {
        if (!processedPhases.found(phase.name()))
        {
            const scalar mass =
                gSum(alphaRho_[phaseIndex(phase.name())].primitiveField()
                   * mesh_.V().field());
            Info << "    " << phase.name() << ": " << mass << " kg" << endl;
        }
    }


    calcAlphas();

    return tPCR;
}


// ************************************************************************* //