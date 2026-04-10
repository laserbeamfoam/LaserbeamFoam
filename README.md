# `laserbeamFoam` solvers

<p align="center">
  <img src="images/LBF.gif" alt="LBF gif" style="width:1000px;">
</p>

## Overview

Presented here is a growing suite of solvers that describe the laser-substrate
 interaction. This repository begins with the `laserbeamFoam` solver. Additional
 solvers are being added incrementally.

Currently, this repository contains two solvers:

### `laserbeamFoam`

A volume-of-fluid (VOF) solver for studying high-energy-density laser-based
advanced manufacturing processes and laser-substrate interactions. This
implementation treats the metallic substrate and shielding gas phase as
in-compressible. The solver fully captures the metallic substrate's
fusion/melting state transition. For the vapourisation of the substrate, the
explicit volumetric dilation due to the vapourisation state transition is
neglected; instead, a phenomenological recoil pressure term is used to capture
the contribution to the momentum and energy fields due to vaporisation events.
laserbeamFoam also captures surface tension effects, the temperature dependence
of surface tension (Marangoni) effects, latent heat effects due to
melting/fusion (and vapourisation), buoyancy effects due to the thermal
expansion of the phases using a Boussinesq approximation, and momentum damping
due to solidification.
A ray-tracing algorithm is implemented that permits the incident Gaussian laser
beam to be discretised into several 'Rays' based on the computational grid
resolution. The 'Rays' of this incident laser beam are then tracked through the
domain through their multiple reflections, with the energy deposited by each
ray determined through the Fresnel equations. The solver approach is extended
from the adiabatic two-phase interFoam code developed by
[OpenCFD Ltd.](http://openfoam.com/) to include non-isothermal state transition
physics and ray-tracing heat source application.

### `compressiblelaserbeamFoam`

An extension of the `laserbeamfoam` solver to multi-component metallic
substrates. This solver can simulate (2N+1)-component systems systems where the 2N is because N-components can exist in the domain in their condensed and vapourised states and the mass transfer and volumetric dilation between these states is fully captured. Diffusion is treated through a Fickian diffusion model
with the diffusivity specified through 'diffusion pairs', and the interface
compression is again specified pair-wise. The miscible phases in the simulation
should have diffusivity specified between them, and immiscible phase pairs
should have an interface compression term specified between them (typically 1).

Target applications for the solvers included in this repository include:

* Laser Welding
* Laser Drilling
* Laser Powder Bed Fusion
* Selective Laser Melting

## Installation

The `OpenFoam_com_main` branch compiles with openfoam v2506, while the `Openfoam_org_main` branch compiles
 with OpenFOAM10. To install the `laserbeamFoam` solvers, first, install and
 load a compatible version of OpenFOAM, then clone and build the `laserbeamFoam`
 library:

```bash
https://github.com/laserbeamfoam/LaserbeamFoam.git
cd LaserbeamFoam && ./Allwmake -j
```

where the `-j` option uses all CPU cores available for building.

The installation can be tested using the tutorial cases described below.

### Optional: Installation of the LIGGGHTS® Discrete Element Model Solver

Some of the tutorial cases use a discrete element method (DEM) solver called
 LIGGGHTS to simulate the creation of a powder bed, e.g. see [this powder bed
 fusion tutorial](tutorials/laserbeamFoam/LPBF_small/README.md).
 For these cases, if available, the `liggghts` executable will be used in the
 case pre-processing process.

On Linux, LIGGGHTS® can be installed with

```bash
# Install required dependencies
sudo apt update
sudo apt install -y build-essential cmake gfortran git \
  libfftw3-dev libjpeg-dev libpng-dev libvtk6-dev \
  libopenmpi-dev openmpi-bin
  
# Clone the LIGGGHTS repository
git clone https://github.com/CFDEMproject/LIGGGHTS-PUBLIC.git
cd LIGGGHTS-PUBLIC/src

# Compile
make auto

# The `liggghts` executable should now be available in this directory
```

While on macOS, LIGGGHTS® can be installed with

```bash
# Install required dependencies using Homebrew
brew install cmake gcc openmpi vtk
  
# Clone the LIGGGHTS repository
git clone https://github.com/CFDEMproject/LIGGGHTS-PUBLIC.git

# Compile
# You may need to update the vtk version in the cmake command to the version
# installed on your system (i.e., replace 9.4.2_1 with another version)
cd LIGGGHTS-PUBLIC
mkdir build
cd build
cmake ../src -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DVTK_DIR=/opt/homebrew/Cellar/vtk/9.4.2_1/lib/cmake/vtk-9.4
make

# The `liggghts` executable should now be available in this directory
```

For convenience, you can add add `liggghts` to your PATH (e.g. in `~/.bashrc`):

```bash
export PATH="~/LIGGGHTS-PUBLIC/build:$PATH"
```

where the location should be updated to match the location on your system.

## Tutorial Cases

The tutorial cases can be run with the included `Allrun` scripts, i.e.

```bash
./Allrun
```

The `Allrun` script prepares the mesh and fields, and runs the solver. Typically
 the following steps are performed:

```bash
# Create the 0 directory
cp -r initial 0

# Create the mesh
blockMesh

# Set the initial fields
setFields

# Run the solver in serial
laserbeamFoam

# Or run the solver in parallel, e.g. on 6 cores
#decomposePar
#mpirun -np 6 laserbeamFoam -parallel &> log.laserbeamFoam
```

Cases can be cleaned and reset using the included `Allclean` scripts, i.e.

```bash
./Allclean
```


## Algorithm

Initially, the solver loads the mesh, reads in fields and boundary conditions,
 and selects the turbulence model (if specified). The main solver loop is then
 initiated. First, the time step is dynamically modified to ensure numerical
 stability. Next, the two-phase fluid mixture properties and turbulence
 quantities are updated. The discretised phase-fraction equation is then solved
 for a user-defined number of subtime steps (typically 3) using the
 multidimensional universal limiter with explicit solution solver [MULES](https://openfoam.org/release/2-3-0/multiphase/).
 This solver is included in the OpenFOAM library and performs conservative
 solutions of hyperbolic convective transport equations with defined bounds (0
 and 1 for $α_1$). Once the updated phase field is obtained, the program enters
 the pressure–velocity loop, where p and u are corrected alternatingly. $T$ is
 also solved in this loop, so that the buoyancy predictions are correct for the
 $U$ and $p$ fields. Correcting the pressure and velocity fields in the sequence
 is known as pressure implicit with the splitting of operators (PISO). In the
 OpenFOAM environment, PISO is repeated for multiple iterations at each time
 step. This process is called the merged PISO- semi-implicit method for
 pressure-linked equations (SIMPLE) or the pressure-velocity loop (PIMPLE)
 process, where SIMPLE is an iterative pressure–velocity solution algorithm.
 PIMPLE continues for a user-specified number of iterations.

The main solver loop iterates until program termination. A summary of the
 simulation algorithm is presented below:

* `laserbeamFoam` Simulation Algorithm Summary:

  * Initialise simulation data and mesh

  * WHILE $t < t_{\text{end}}$ DO

    1. Update $\Delta t$ for stability

    2. Phase equation sub-cycle

    3. Update interface location for the heat source application

    4. Update fluid properties

    5. Ray-tracing for Heat Source application at the surface

    6. PISO Loop

        1. Form $U$ equation

        2. Energy Transport Loop

            1. Solve $T$ equation  
            2. Update fluid fraction field  
            3. Re-evaluate source terms due to latent heat

        3. PISO

            1. Obtain and correct face fluxes  
            2. Solve $p$ Poisson equation  
            3. Correct $U$

    7. Write fields

There are no constraints on how the computational domain is discretised.

## Visualising the rays in ParaView

`laserbeamFoam` writes the individual ray beams to `VTK/rays_<LASER_NAME>_<TIME>.vtk`,
 where `<LASER_NAME>` is the laser's name within the `laser` sub-dict of `constant/LaserProperties`
 and `<TIME>` is the time value. ParaView recognises that these files are in a
 sequence, so they can all be loaded together: `File` -> `Open...` -> Select
 `rays_..vtk`. As the VTK files do not store time-step information, by default,
 ParaView assumes the time-step size for the rays is 1 s; however, you can use
 the ParaView “Temporal Shift Scale” filter on the rays object to sync the ray
 time with the OpenFOAM model time, where the OpenFOAM time-step value (e.g.
 1e-5) is used as the `Scale`. For convenience, a `VTK/rays_<LASER_NAME>.vtk.series`
 files is written at the end of the simulation. This `.vtk.series` file can be
 opened directly in ParaView, which loads the VTK files at the correct physical
 times (no need for the “Temporal Shift Scale” filter). Once the ray VTKs are
 loaded, they can be viewed by setting their colour appropriately (e.g. to
 black); in addition, it is often convenient to increase the "Line Width" or
 use the "Tube" filter.

## License

OpenFOAM, and by extension, the `laserbeamFoam` application, is licensed free
 and open source only under the [GNU General Public Licence version 3](https://www.gnu.org/licenses/gpl-3.0.en.html).
 One reason for OpenFOAM’s popularity is that its users are granted the freedom
 to modify and redistribute the software and have a right to continued free use
 within the terms of the GPL.

## Acknowledgements

Tom Flint and Joe Robson thank the EPSRC for financial support through the
 associated programme grant LightFORM (EP/R001715/1). Joe Robson thanks the
 Royal Academy of Engineering/DSTL for funding through the RAEng/DSTL Chair in
 Alloys for Extreme Environments.

Philip Cardiff and Gowthaman Parivendhan authors gratefully acknowledge financial
 support from I-Form, funded by Science Foundation Ireland (SFI) Grant Numbers
 16/RC/3872 and 21/RC/10295 P2, co-funded under the European Regional Development
 Fund and by I-Form industry partners. In addition, Philip Cardiff received
 funding from the European Research Council (ERC) under the European Union’s
 Horizon 2020 research and innovation programme (Grant Agreement No. 101088740),
 and acknowledges financial support from the Irish Research Council through the
 Laureate programme, grant number IRCLA/2017/45, and Bekaert, through the Bekaert
 University Technology Centre (UTC) at University College Dublin
 [www.ucd.ie/bekaert](www.ucd.ie/bekaert).

## Citing This Work

If you use `laserbeamFoam` in your work. Please use the following to cite our work:

```bibtex
laserbeamFoam: Laser ray-tracing and thermally induced state transition
simulation toolkit. TF Flint, JD Robson, G Parivendhan, P Cardiff - SoftwareX,
2023 - https://doi.org/10.1016/j.softx.2022.101299
```

## References

Flint, T. F., Robson, J. D., Parivendhan, G., & Cardiff, P. (2023).
 laserbeamFoam: Laser ray-tracing and thermally induced state transition
 simulation toolkit. SoftwareX, 21, 101299.

Flint, T. F., Parivendhan, G., Ivankovic, A., Smith, M. C., & Cardiff, P. (2022).
 beamWeldFoam: Numerical simulation of high energy density fusion and
 vapourisation-inducing processes. SoftwareX, 18, 101065.

Flint, T. F., et al. A fundamental analysis of factors affecting chemical
 homogeneity in the laser powder bed fusion process. International Journal of
 Heat and Mass Transfer 194 (2022): 122985.

Flint, T. F., T. Dutilleul, and W. Kyffin. A fundamental investigation into the
 role of beam focal point, and beam divergence, on thermo-capillary stability
 and evolution in electron beam welding applications. International Journal of
 Heat and Mass Transfer 212 (2023): 124262.

Parivendhan, G., Cardiff, P., Flint, T., Tuković, Ž., Obeidi, M., Brabazon, D.,
 Ivanković, A. (2023) A numerical study of processing parameters and their
 effect on the melt-track profile in Laser Powder Bed Fusion processes, Additive
 Manufacturing, 67, 10.1016/j.addma.2023.103482.

## Disclaimer

This offering is not approved or endorsed by OpenCFD Limited, producer and
 distributor of the OpenFOAM software via [www.openfoam.com](https://www.openfoam.com),
 and owner of the OPENFOAM® and OpenCFD® trade marks.

## Acknowledgement

OPENFOAM® is a registered trademark of OpenCFD Limited, producer and distributor
 of the OpenFOAM software via [www.openfoam.com](https://www.openfoam.com).

![visitors](https://visitor-badge.deta.dev/badge?page_id=micmog.LaserbeamFoam)
