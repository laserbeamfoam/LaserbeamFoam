#ifndef solidElasticity_H
#define solidElasticity_H

#include "fvCFD.H"

namespace Foam
{

class solidElasticity
{
    const fvMesh& mesh_;
    const volScalarField& alphaMetal_;
    const volScalarField& epsilon1_;
    const volScalarField& T_;

    volScalarField& fSolid_;
    volVectorField& D_;
    volSymmTensorField& epsilonElast_;
    volSymmTensorField& sigmaElast_;

    dimensionedScalar E_;
    dimensionedScalar nu_;
    dimensionedScalar alphaT_;
    dimensionedScalar Tref_;
    dimensionedScalar solidDamping_;
    dimensionedScalar fSolidMin_;

    label maxIter_;
    scalar tolerance_;
    scalar relTol_;
    word coupling_;
    bool compactNormalStress_;

    void readDict(const dictionary& solidDict);

public:
    solidElasticity
    (
        const fvMesh& mesh,
        const volScalarField& alphaMetal,
        const volScalarField& epsilon1,
        const volScalarField& T,
        volScalarField& fSolid,
        volVectorField& D,
        volSymmTensorField& epsilonElast,
        volSymmTensorField& sigmaElast,
        const dictionary& solidDict
    );

    void read(const dictionary& solidDict);

    void updateFSolid();

    void solve();

    const dimensionedScalar& solidDamping() const;
};

} // End namespace Foam

#endif
