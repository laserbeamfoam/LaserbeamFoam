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

#include "temperatureTable.H"
#include "messageStream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::temperatureTable::readControls(const dictionary& controls)
{
    // defualts cover anything a laser case can reach
    Tmin_ = controls.getOrDefault<scalar>("Tmin", 0.0);
    Tmax_ = controls.getOrDefault<scalar>("Tmax", 1.0e4);
    deltaT_ = controls.getOrDefault<scalar>("deltaT", 1.0);

    if (!(Tmax_ > Tmin_) || !(deltaT_ > 0))
    {
        FatalIOErrorInFunction(controls)
            << "propertyTables: need Tmax > Tmin and deltaT > 0, got Tmin = "
            << Tmin_ << ", Tmax = " << Tmax_ << ", deltaT = " << deltaT_
            << exit(FatalIOError);
    }

    rDeltaT_ = 1.0/deltaT_;
}


Foam::label Foam::temperatureTable::nSamples() const
{
    const label n = label(round((Tmax_ - Tmin_)/deltaT_)) + 1;

    if (n < 2 || n > 20000000)
    {
        FatalErrorInFunction
            << "propertyTables: (Tmax - Tmin)/deltaT + 1 = " << n
            << " samples for " << name_ << "; expected between 2 and 2e7"
            << exit(FatalError);
    }

    return n;
}


void Foam::temperatureTable::sample(const Function1<scalar>& f)
{
    values_.setSize(nSamples());

    forAll(values_, i)
    {
        values_[i] = f.value(Tmin_ + i*deltaT_);
    }
}


void Foam::temperatureTable::sample(const Polynomial<8>& p)
{
    values_.setSize(nSamples());

    forAll(values_, i)
    {
        values_[i] = p.value(Tmin_ + i*deltaT_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperatureTable::temperatureTable
(
    const word& keyword,
    const dictionary& dict,
    const dictionary& controls,
    const word& legacyPolynomialKeyword
)
:
    name_(keyword),
    description_(),
    Tmin_(0),
    Tmax_(0),
    deltaT_(1),
    rDeltaT_(1),
    values_()
{
    readControls(controls);

    if (dict.found(keyword))
    {
        if
        (
            legacyPolynomialKeyword != word::null
         && dict.found(legacyPolynomialKeyword)
        )
        {
            Info<< "temperatureTable: " << dict.dictName() << "/" << keyword
                << " given as a Function1 and as " << legacyPolynomialKeyword
                << ": using " << keyword << endl;
        }

        autoPtr<Function1<scalar>> f(Function1<scalar>::New(keyword, dict));

        description_ = string("Function1 ") + (*f).type();

        sample(*f);
    }
    else if
    (
        legacyPolynomialKeyword != word::null
     && dict.found(legacyPolynomialKeyword)
    )
    {
        // old style a0 + a1 T + a2 T^2 ... coeficients
        const Polynomial<8> p(dict.lookup(legacyPolynomialKeyword));

        description_ = string("polynomial ") + legacyPolynomialKeyword;

        sample(p);
    }
    else
    {
        string alternatives
        (
            "as a Function1 of the temperature "
            "(constant, table, tableFile, polynomial, ...)"
        );

        if (legacyPolynomialKeyword != word::null)
        {
            alternatives +=
                " or as the polynomial coefficients " + legacyPolynomialKeyword;
        }

        FatalIOErrorInFunction(dict)
            << "Property " << keyword << " not found in dictionary "
            << dict.dictName() << ": give it " << alternatives.c_str()
            << exit(FatalIOError);
    }
}


Foam::temperatureTable::temperatureTable
(
    const word& name,
    const scalar value,
    const dictionary& controls
)
:
    name_(name),
    description_("constant"),
    Tmin_(0),
    Tmax_(0),
    deltaT_(1),
    rDeltaT_(1),
    values_()
{
    readControls(controls);

    values_.setSize(nSamples(), value);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::temperatureTable::isConstant() const
{
    return min(values_) == max(values_);
}


void Foam::temperatureTable::values
(
    const UList<scalar>& T,
    UList<scalar>& result
) const
{
    if (result.size() != T.size())
    {
        FatalErrorInFunction
            << "Size mismatch for " << name_ << ": " << T.size()
            << " temperatures, " << result.size() << " results"
            << abort(FatalError);
    }

    forAll(T, i)
    {
        result[i] = value(T[i]);
    }
}


Foam::tmp<Foam::scalarField> Foam::temperatureTable::values
(
    const UList<scalar>& T
) const
{
    auto tresult = tmp<scalarField>::New(T.size());

    values(T, tresult.ref());

    return tresult;
}


void Foam::temperatureTable::printInfo(const UList<scalar>& Tsample) const
{
    Info<< "    " << name_ << ": " << description_.c_str()
        << (isConstant() ? " (constant)" : "")
        << ", cached " << Tmin_ << ".." << Tmax_ << " K every " << deltaT_
        << " K";

    if (!isConstant())
    {
        Info<< ", values:";

        forAll(Tsample, i)
        {
            Info<< " " << Tsample[i] << "K:" << value(Tsample[i]);
        }
    }
    else
    {
        Info<< ", value " << values_[0];
    }

    Info<< endl;
}


// ************************************************************************* //
