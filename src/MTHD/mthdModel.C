/*---------------------------------------------------------------------------*\
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
    along with laserbeamFoam.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "mthdModel.H"

#include "fvCFD.H"
#include "zeroGradientFvPatchFields.H"

namespace Foam
{

defineTypeNameAndDebug(mthdModel, 0);

dimensionedScalar mthdModel::lookupElectricalConductivity
(
    const dictionary& targetDictionary
)
{
    static const dimensionSet conductivityDimensions(-1, -3, 3, 0, 0);
    static const word oldKeyword("elec_resistivity");
    static const word newKeyword("Elec_conductivity");

    if
    (
        targetDictionary.found(newKeyword)
     && targetDictionary.found(oldKeyword)
    )
    {
        const dimensionedScalar conductivity
        (
            newKeyword,
            conductivityDimensions,
            targetDictionary
        );

        const scalar oldScalarValue = targetDictionary.get<scalar>(oldKeyword);
        const scalar conductivityFromOld = 1.0/oldScalarValue;
        constexpr scalar toleranceForCheck = 10.0 * Foam::SMALL;

        Info<< "Warning: using deprecated keyword `" << oldKeyword << "`."
            << " Please use `" << newKeyword << "` in future." << nl
            << "Warning: detected both " << oldKeyword << " and "
            << newKeyword << " in " << targetDictionary.name() << nl
            << "Prioritising " << newKeyword << "." << endl;

        if (mag(conductivity.value() - conductivityFromOld) > toleranceForCheck)
        {
            FatalErrorInFunction
                << "CAUTION: " << oldKeyword << " and " << newKeyword
                << " values are dissimilar in " << targetDictionary.name()
                << nl
                << "Tolerance   = " << toleranceForCheck << nl
                << oldKeyword << "^-1 = " << conductivityFromOld << nl
                << newKeyword << " = " << conductivity.value() << nl
                << exit(FatalError);
        }

        return conductivity;
    }

    if (targetDictionary.found(oldKeyword))
    {
        Info<< "Warning: using deprecated keyword `" << oldKeyword << "`."
            << " Please use `" << newKeyword << "` in future." << endl;

        return dimensionedScalar
        (
            newKeyword,
            conductivityDimensions,
            1.0/targetDictionary.get<scalar>(oldKeyword)
        );
    }

    return dimensionedScalar
    (
        newKeyword,
        conductivityDimensions,
        targetDictionary
    );
}


dimensionedScalar mthdModel::readElectricalConductivity
(
    const dictionary& targetDictionary
)
{
    return lookupElectricalConductivity(targetDictionary);
}


Switch mthdModel::electromagneticsEnabled
(
    const fvMesh& mesh,
    const word& propertiesDictName
)
{
    return
        IOdictionary
        (
            IOobject
            (
                propertiesDictName,
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        ).lookupOrDefault<Switch>("electromagnetics", false);
}


mthdModel::mthdModel(fvMesh& mesh, const word& propertiesDictName)
:
    mesh_(mesh),
    bpiso_(mesh_, "BPISO"),
    phiH_
    (
        IOobject
        (
            "phiH",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ
        ).typeHeaderOk<surfaceScalarField>(true)
      ? surfaceScalarField
        (
            IOobject
            (
                "phiH",
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_
        )
      : surfaceScalarField
        (
            IOobject
            (
                "phiH",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::flux
            (
                volVectorField
                (
                    IOobject
                    (
                        "H",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_
                )
            )
        )
    ),
    H_
    (
        IOobject
        (
            "H",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    pH_
    (
        IOobject
        (
            "pH",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    B_
    (
        IOobject
        (
            "B",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("B", dimensionSet(1, 0, -2, 0, 0), vector::zero)
    ),
    J_
    (
        IOobject
        (
            "J",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("J", dimensionSet(0, -2, 0, 0, 0), vector::zero)
    ),
    divB_
    (
        IOobject
        (
            "divB",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("divB", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    eleccond_
    (
        IOobject
        (
            "eleccond",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("eleccond", dimensionSet(-1, -3, 3, 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    magperm_
    (
        IOobject
        (
            "magperm",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("magperm", dimensionSet(1, 1, -2, 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    oneoversigma_
    (
        IOobject
        (
            "oneoversigma",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("oneoversigma", dimensionSet(1, 3, -3, 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Joulesource_
    (
        IOobject
        (
            "Joulesource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Joulesource", dimensionSet(1, -1, -3, 0, 0), 0.0)
    ),
    JcrossB_
    (
        IOobject
        (
            "JcrossB",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("JcrossB", dimensionSet(1, -2, -2, 0, 0), vector::zero),
        zeroGradientFvPatchVectorField::typeName
    )
{
    (void)propertiesDictName;

    if
    (
        IOobject
        (
            "phiH",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ
        ).typeHeaderOk<surfaceScalarField>(true)
    )
    {
        Info<< "Reading face flux ";
    }
    else
    {
        Info<< "Calculating face flux ";
    }

    Info<< phiH_.name() << nl << endl;

    H_.correctBoundaryConditions();
    pH_.correctBoundaryConditions();
    B_.write();
}


void mthdModel::correctIncompressibleProperties
(
    const volScalarField& alpha1,
    const dictionary& transportPropertiesMetal,
    const dictionary& transportPropertiesGas
)
{
    const dimensionedScalar eleccond1 =
        lookupElectricalConductivity(transportPropertiesMetal);
    const dimensionedScalar eleccond2 =
        lookupElectricalConductivity(transportPropertiesGas);

    const dimensionedScalar magperm1
    (
        "Mag_permitivity",
        dimensionSet(1, 1, -2, 0, 0),
        transportPropertiesMetal
    );

    const dimensionedScalar magperm2
    (
        "Mag_permitivity",
        dimensionSet(1, 1, -2, 0, 0),
        transportPropertiesGas
    );

    eleccond_ = alpha1*eleccond1 + (1.0 - alpha1)*eleccond2;
    magperm_ = alpha1*magperm1 + (1.0 - alpha1)*magperm2;

    eleccond_.correctBoundaryConditions();
    magperm_.correctBoundaryConditions();
}


void mthdModel::resetTransportProperties()
{
    eleccond_ *= 0.0;
    magperm_ *= 0.0;
}


void mthdModel::accumulateTransportProperties
(
    const volScalarField& alpha,
    const dimensionedScalar& electricalConductivity,
    const dimensionedScalar& magneticPermeability
)
{
    eleccond_ += alpha*electricalConductivity;
    magperm_ += alpha*magneticPermeability;
}


void mthdModel::finaliseTransportProperties()
{
    eleccond_.correctBoundaryConditions();
    magperm_.correctBoundaryConditions();
}


void mthdModel::updateLorentzForce()
{
    JcrossB_ = J_ ^ (magperm_*H_);
    JcrossB_.correctBoundaryConditions();
}


void mthdModel::solve
(
    const surfaceScalarField& phi,
    const volVectorField& U
)
{
    while (bpiso_.correct())
    {
        oneoversigma_ = 1.0/eleccond_;
        oneoversigma_.correctBoundaryConditions();

        surfaceScalarField muPhiH(fvc::interpolate(magperm_)*phiH_);
        surfaceScalarField muPhi(fvc::interpolate(magperm_)*phi);

        volTensorField TgradH(fvc::grad(H_)().T());
        TgradH.correctBoundaryConditions();

        fvVectorMatrix HEqn
        (
            fvm::ddt(magperm_, H_)
          + fvm::div(muPhi, H_)
          - fvc::div(muPhiH, U)
          - fvm::laplacian(oneoversigma_, H_)
          + fvc::div(oneoversigma_*TgradH)
        );

        HEqn.relax();
        HEqn.solve();

        H_.correctBoundaryConditions();

        volScalarField rAB(1.0/HEqn.A());
        rAB.correctBoundaryConditions();
        surfaceScalarField rABf("rABf", fvc::interpolate(rAB));

        volVectorField HbyAH("HbyAH", H_);
        HbyAH = rAB*HEqn.H();

        surfaceScalarField phiHbyAH("phiHbyAH", fvc::flux(HbyAH));

        while (bpiso_.correctNonOrthogonal())
        {
            fvScalarMatrix pHEqn
            (
                fvm::laplacian(rABf, pH_)
             == fvc::div(fvc::interpolate(magperm_)*phiHbyAH)
            );

            pHEqn.solve();
            pH_.correctBoundaryConditions();

            if (bpiso_.finalNonOrthogonalIter())
            {
                phiH_ = phiHbyAH - pHEqn.flux()/fvc::interpolate(magperm_);
                H_.correctBoundaryConditions();
            }
        }

        Info<< "H magnetic flux divergence error {div(mu H)}= "
            << mesh_.time().deltaTValue()
               *mag(fvc::div(fvc::interpolate(magperm_)*phiH_))()
                    .weightedAverage(mesh_.V()).value()
            << endl;

        B_ = magperm_*H_;
        J_ = fvc::curl(H_);
        J_.correctBoundaryConditions();

        Joulesource_ = (1.0/eleccond_)*magSqr(J_);
        Joulesource_.correctBoundaryConditions();

        divB_ = fvc::div(fvc::interpolate(magperm_)*phiH_);
        divB_.correctBoundaryConditions();
    }
}

} // End namespace Foam
