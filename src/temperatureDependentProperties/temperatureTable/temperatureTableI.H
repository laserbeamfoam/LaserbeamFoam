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

inline Foam::scalar Foam::temperatureTable::value(const scalar T) const
{
    const scalar x = (T - Tmin_)*rDeltaT_;

    // below the table (also catches a NaN temperture)
    if (!(x > 0))
    {
        return values_[0];
    }

    const label last = values_.size() - 1;

    if (x >= last)
    {
        return values_[last];
    }

    const label i = static_cast<label>(x);
    const scalar f = x - i;

    return values_[i] + f*(values_[i + 1] - values_[i]);
}


// ************************************************************************* //
