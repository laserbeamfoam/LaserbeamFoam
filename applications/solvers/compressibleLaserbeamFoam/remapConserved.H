// After an AMR change: the mapper volume-averages c_i and rhoM (mass kept), and T, U and epsilon1 are
// recovered from the mapped energy, momentum and latent heat, so merging cells conserves those too.
// mesh.update() runs before anything is solved this step, so the old-time values are reset as well
{
    // rhoCv and CvMass from the mapped c_i and rhoM, exactly as TEqn.H builds them
    #include "updateRhoCv.H"

    T = rhoCvTAmr/rhoCv;
    T.correctBoundaryConditions();
    T.oldTime() == T;

    U = rhoMUAmr/rhoM;
    U.correctBoundaryConditions();
    U.oldTime() == U;

    // mass-weighted, so the gas cells' epsilon1 = 1 does not melt a merged solid cell
    epsilon1 = rhoMEps1Amr/rhoM;
    epsilon1.correctBoundaryConditions();
    epsilon1.oldTime() == epsilon1;

    // TEqn.H rebuilds rhoCv.oldTime() = rhoM.oldTime()*CvMass.oldTime() and the main
    // loop sets rhoM.oldTime() = rhoM, so the step starts from the mapped energy
    CvMass.oldTime() == CvMass;
}
