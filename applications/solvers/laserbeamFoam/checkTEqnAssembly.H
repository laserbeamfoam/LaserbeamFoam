// Whether the transient, convection and conduction terms of the energy
// equation and their continuity correction can be assembled once per time
// step rather than in every correction of the liquid fraction (TEqn.H): their
// coefficients do not depend on the temperature of the current correction
// when the convection scheme is upwind (the weights depend only on the flux)
// and the conduction term has no non-orthogonal correction (an uncorrected or
// orthogonal scheme, or a mesh whose non-orthogonal correction vectors are
// zero). Checked in every time step, as the schemes may be re-read; the
// non-orthogonal correction vectors only at start-up and when the mesh
// changes.
{
    if (!TEqnAssemblyChecked || mesh.changing())
    {
        TEqnNonOrthCorrectionZero =
            max(mag(mesh.nonOrthCorrectionVectors())).value() == 0;
    }

    ITstream& divStream = mesh.divScheme("div(rhophicp,T)");
    const word divDiscretisation(divStream);
    const word divInterpolation(divStream);

    ITstream& lapStream = mesh.laplacianScheme("laplacian(kappa,T)");
    const word lapDiscretisation(lapStream);
    const word lapInterpolation(lapStream);
    const word lapSnGrad(lapStream);

    const bool assembleOnce =
        divDiscretisation == "Gauss" && divInterpolation == "upwind"
     && (
            lapSnGrad == "uncorrected"
         || lapSnGrad == "orthogonal"
         || TEqnNonOrthCorrectionZero
        );

    if (assembleOnce != TEqnAssembleOnce || !TEqnAssemblyChecked)
    {
        Info<< "TEqn: the temperature-independent terms are assembled "
            << (assembleOnce ? "once per time step" : "in every correction")
            << endl;
    }

    TEqnAssembleOnce = assembleOnce;
    TEqnAssemblyChecked = true;
}
