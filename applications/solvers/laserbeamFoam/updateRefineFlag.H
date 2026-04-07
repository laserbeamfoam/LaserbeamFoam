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
    const volScalarField& alphaField
)
{
    const scalar gradCondThresh = 0.5;
    const scalar liquidStateMin = -0.5;
    const scalar liquidStateMax = 1.5;

    refineFlag *= 0.0;

    scalarField& refineInternal = refineFlag.ref();
    const scalarField& stateInternal = State.internalField();
    const tmp<volScalarField> tGradCondMag = mag(fvc::grad(alphaField));
    const scalarField& gradCondMagInternal = tGradCondMag().internalField();

    forAll(stateInternal, celli)
    {
        // Refine liquid and mushy condensed material:
        //   State = 0 -> liquid
        //   State = 1 -> mushy
        const bool liquidVolume =
            stateInternal[celli] > liquidStateMin
         && stateInternal[celli] < liquidStateMax;

        // Refine the gas/condensed interface using the condensate gradient.
        const bool condensedInterface =
            gradCondMagInternal[celli] > gradCondThresh;

        if
        (
            liquidVolume
         || condensedInterface
        )
        {
            refineInternal[celli] = 1.0;
        }
    }

    refineFlag.correctBoundaryConditions();
}

} // namespace Foam

#endif
