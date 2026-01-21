#include "solidElasticity.H"

namespace Foam
{

solidElasticity::solidElasticity
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
)
:
    mesh_(mesh),
    alphaMetal_(alphaMetal),
    epsilon1_(epsilon1),
    T_(T),
    fSolid_(fSolid),
    D_(D),
    epsilonElast_(epsilonElast),
    sigmaElast_(sigmaElast),
    E_("E", dimensionSet(1, -1, -2, 0, 0), 2.0e11),
    nu_("nu", dimless, 0.30),
    alphaT_("alphaT", dimensionSet(0, 0, 0, -1, 0), 1.2e-5),
    Tref_("Tref", dimensionSet(0, 0, 0, 1, 0), 300.0),
    solidDamping_("solidDamping", dimensionSet(1, -3, -1, 0, 0), 1.0e14),
    fSolidMin_("fSolidMin", dimless, 1e-6),
    maxIter_(10),
    tolerance_(1e-6),
    relTol_(0),
    coupling_("implicit")
{
    readDict(solidDict);
    updateFSolid();
}

void solidElasticity::readDict(const dictionary& solidDict)
{
    if (!solidDict.found("solidElasticity"))
    {
        return;
    }

    const dictionary& elastDict = solidDict.subDict("solidElasticity");

    coupling_ = elastDict.lookupOrDefault<word>("coupling", "implicit");
    if (coupling_ != "implicit")
    {
        FatalErrorInFunction
            << "Only implicit coupling is supported. Found coupling = "
            << coupling_ << nl << exit(FatalError);
    }

    maxIter_ = elastDict.lookupOrDefault<label>("maxIter", 10);
    tolerance_ = elastDict.lookupOrDefault<scalar>("tolerance", 1e-6);
    relTol_ = elastDict.lookupOrDefault<scalar>("relTol", 0);

    E_ = dimensionedScalar
    (
        "E",
        dimensionSet(1, -1, -2, 0, 0),
        elastDict.lookupOrDefault<scalar>("E", 2.0e11)
    );

    nu_ = dimensionedScalar
    (
        "nu",
        dimless,
        elastDict.lookupOrDefault<scalar>("nu", 0.30)
    );

    alphaT_ = dimensionedScalar
    (
        "alphaT",
        dimensionSet(0, 0, 0, -1, 0),
        elastDict.lookupOrDefault<scalar>("alphaT", 1.2e-5)
    );

    Tref_ = dimensionedScalar
    (
        "Tref",
        dimensionSet(0, 0, 0, 1, 0),
        elastDict.lookupOrDefault<scalar>("Tref", 300.0)
    );

    solidDamping_ = dimensionedScalar
    (
        "solidDamping",
        dimensionSet(1, -3, -1, 0, 0),
        elastDict.lookupOrDefault<scalar>("solidDamping", 1.0e14)
    );

    fSolidMin_ = dimensionedScalar
    (
        "fSolidMin",
        dimless,
        elastDict.lookupOrDefault<scalar>("fSolidMin", 1e-6)
    );
}

void solidElasticity::read(const dictionary& solidDict)
{
    if (solidDict.found("solidElasticity"))
    {
        readDict(solidDict);
    }
}

void solidElasticity::updateFSolid()
{
    fSolid_ = max(min(alphaMetal_, scalar(1.0)), scalar(0.0))
        *(scalar(1.0) - max(min(epsilon1_, scalar(1.0)), scalar(0.0)));

    fSolid_.correctBoundaryConditions();
}

void solidElasticity::solve()
{
    updateFSolid();

    const dimensionedScalar mu0
    (
        "mu0",
        dimensionSet(1, -1, -2, 0, 0),
        E_.value()/(2.0*(1.0 + nu_.value()))
    );

    const dimensionedScalar lambda0
    (
        "lambda0",
        dimensionSet(1, -1, -2, 0, 0),
        E_.value()*nu_.value()/((1.0 + nu_.value())*(1.0 - 2.0*nu_.value()))
    );

    const dimensionedScalar threeK0
    (
        "threeK0",
        dimensionSet(1, -1, -2, 0, 0),
        E_.value()/(1.0 - 2.0*nu_.value())
    );

    const volScalarField fSolidEff
    (
        IOobject
        (
            "fSolidEff",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        max(fSolid_, fSolidMin_)
    );

    const volScalarField mu
    (
        IOobject
        (
            "muElast",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fSolidEff*mu0
    );

    const volScalarField lambda
    (
        IOobject
        (
            "lambdaElast",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fSolidEff*lambda0
    );

    const volScalarField threeKalpha
    (
        IOobject
        (
            "threeKalpha",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fSolidEff*threeK0*alphaT_
    );

    scalar initialResidual = GREAT;
    scalar initialResidual0 = GREAT;
    label iter = 0;

    do
    {
        fvVectorMatrix DEqn
        (
            fvm::laplacian(2*mu + lambda, D_, "laplacian(DD,D)")
         ==
            fvc::grad(threeKalpha*(T_ - Tref_))
        );

        if (mesh_.nCells() > 0)
        {
            DEqn.setComponentReference(0, 0, vector::X, 0);
            DEqn.setComponentReference(0, 0, vector::Y, 0);
            DEqn.setComponentReference(0, 0, vector::Z, 0);
        }

        DEqn.relax();
        initialResidual = DEqn.solve().max().initialResidual();

        if (iter == 0)
        {
            initialResidual0 = max(initialResidual, VSMALL);
        }

        ++iter;
    }
    while
    (
        iter < maxIter_
     && initialResidual > tolerance_
     && (relTol_ <= 0 || initialResidual/initialResidual0 > relTol_)
    );

    volTensorField gradD(fvc::grad(D_));

    epsilonElast_ = symm(gradD) - (alphaT_*(T_ - Tref_))*I;
    sigmaElast_ = mu*twoSymm(gradD) + (lambda*I)*tr(gradD)
        - (threeKalpha*(T_ - Tref_))*I;

    epsilonElast_.correctBoundaryConditions();
    sigmaElast_.correctBoundaryConditions();
}

const dimensionedScalar& solidElasticity::solidDamping() const
{
    return solidDamping_;
}

} // End namespace Foam
