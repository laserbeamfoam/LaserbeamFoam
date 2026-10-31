/*---------------------------------------------------------------------------*\
License
    This file is part of laserbeamFoam.

    laserbeamFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    laserbeamFoam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with laserbeamFoam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "temperatureTableViscosity.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(temperatureTableViscosity, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        temperatureTableViscosity,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

// propertyTables controls from the coeffs, else from the parent dict
static dictionary tableControls
(
    const dictionary& viscosityProperties,
    const dictionary& coeffs
)
{
    if (coeffs.found("propertyTables"))
    {
        return coeffs.subDict("propertyTables");
    }

    return viscosityProperties.parent().subOrEmptyDict("propertyTables");
}

}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::viscosityModels::temperatureTableViscosity::calcNu()
{
    const fvMesh& mesh = U_.mesh();

    const volScalarField* Tptr = mesh.findObject<volScalarField>(TName_);

    if (!Tptr)
    {
        WarningInFunction
            << "Temperature field " << TName_ << " not found: evaluating "
            << nu_.name() << " at Tinit = " << Tinit_ << " K" << endl;

        nu_ == dimensionedScalar(dimViscosity, nuTable_.value(Tinit_));

        return;
    }

    const volScalarField& T = *Tptr;

    nuTable_.values(T.primitiveField(), nu_.primitiveFieldRef());

    // boundry values follow the patch temperature
    volScalarField::Boundary& nuBf = nu_.boundaryFieldRef();
    const volScalarField::Boundary& TBf = T.boundaryField();

    forAll(nuBf, patchi)
    {
        nuTable_.values(TBf[patchi], nuBf[patchi]);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::temperatureTableViscosity::temperatureTableViscosity
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    coeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    TName_(coeffs_.getOrDefault<word>("field", "T")),
    Tinit_(coeffs_.getOrDefault<scalar>("Tinit", 300.0)),
    nuTable_("nu", coeffs_, tableControls(viscosityProperties, coeffs_)),
    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh(),
        dimensionedScalar(dimViscosity, Zero)
    )
{
    calcNu();

    Info<< "Selecting temperature dependent viscosity for " << name
        << " from field " << TName_ << nl;
    nuTable_.printInfo(scalarList({300.0, 1000.0, 2000.0, 3000.0, 4000.0}));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::viscosityModels::temperatureTableViscosity::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    coeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");
    TName_ = coeffs_.getOrDefault<word>("field", "T");
    Tinit_ = coeffs_.getOrDefault<scalar>("Tinit", 300.0);

    nuTable_ = Foam::temperatureTable
    (
        "nu",
        coeffs_,
        tableControls(viscosityProperties, coeffs_)
    );

    calcNu();

    return true;
}


// ************************************************************************* //
