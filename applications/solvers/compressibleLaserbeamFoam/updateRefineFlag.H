#ifndef updateRefineFlag_H
#define updateRefineFlag_H

#include "volFields.H"
#include "fvc.H"

namespace Foam
{

inline void updateRefineFlag
(
    volScalarField& refineFlag,
    const volScalarField& State,
    const volScalarField& condensateField
)
{
    const scalar gradCondThresh = 0.5;

    refineFlag *= 0.0;

    scalarField& refineInternal = refineFlag.ref();
    const scalarField& stateInternal = State.internalField();
    const tmp<volScalarField> tGradCondMag = mag(fvc::grad(condensateField));
    const scalarField& gradCondMagInternal = tGradCondMag().internalField();

    forAll(stateInternal, celli)
    {
        if
        (
            stateInternal[celli] > 1.5
         || gradCondMagInternal[celli] > gradCondThresh
        )
        {
            refineInternal[celli] = 1.0;
        }
    }

    refineFlag.correctBoundaryConditions();
}

} // namespace Foam

#endif
