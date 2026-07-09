#ifndef updateRefineFlag_H
#define updateRefineFlag_H

#include "volFields.H"

namespace Foam
{

inline void updateRefineFlag
(
    volScalarField& refineFlag,
    const volScalarField& State,
    const volScalarField& alphaSmoothed
)
{
    const scalar liquidStateMin = -0.5;
    const scalar liquidStateMax = 1.5;
    const scalar interfaceProductThresh = 0.01;

    refineFlag *= 0.0;

    scalarField& refineInternal = refineFlag.ref();
    const scalarField& stateInternal = State.internalField();
    const scalarField& alphaSmoothedInternal = alphaSmoothed.internalField();

    forAll(stateInternal, celli)
    {
        const bool liquidVolume =
            stateInternal[celli] > liquidStateMin
         && stateInternal[celli] < liquidStateMax;

        const bool condensedInterface =
            alphaSmoothedInternal[celli]
          * (1.0 - alphaSmoothedInternal[celli])
            > interfaceProductThresh;

        if (liquidVolume || condensedInterface)
        {
            refineInternal[celli] = 1.0;
        }
    }

    refineFlag.correctBoundaryConditions();
}

} // namespace Foam

#endif