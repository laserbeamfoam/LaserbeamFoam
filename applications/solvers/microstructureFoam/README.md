# microstructureFoam


## Overview

## Installation

The current version of the code utilises the [OpenFoam10 libraries](https://openfoam.org/version/10/). The code has been developed and tested using an Ubuntu installation, but should work on any operating system capable of installing OpenFoam. To install the microstructureFoam solver, first follow the instructions on this page: [OpenFoam 10 Install](https://openfoam.org/download/10-ubuntu/) to install the OpenFoam 10 libraries.



Then navigate to a working folder in a shell terminal, clone the git code repository, and build.

```
$ git clone https://github.com/micmog/laserbeamFoam.git laserbeamFoam
$ cd solver
$ wclean
$ wmake
```
The installation can be tested using the tutorial cases described below.

## Tutorial cases
To run any of the tutorials in serial mode:
```
delete any old simulation files, e.g:
$ rm -r 0* 1* 2* 3* 4* 5* 6* 7* 8* 9*
Then:
$ cp -r initial 0
$ blockMesh
$ setFields
$ TesselateFoam (Optional)
$ microstructureFoam
```
For parallel deployment, using MPI, following the setFields command:
```
$ decomposePar
$ mpirun -np 6 microstructureFoam -parallel >log &
```
for deployment on 6 cores.

## Description of case files
### nucleationDict
| Parameter Name | Value Type | Parameter Description | Unit |
| ------ | ------ | ------ | ------ | 
| `Tu_mean` | dimensionedScalar | Activation undercooling distribution: mean | K |
| `Tu_stdev` | dimensionedScalar | Activation undercooling distribution: standard deviation | K |
| `n_max` | dimensionedScalar| Activation undercooling distribution: site density | [sites/m3] |
| `maxItersNucSet` | scalar | Maximum iterations allowed when setting nucleation sites | - |
| `nucDistFactor` | scalar | Number of interface distances between nucleation sites | - |
| `stopToCheckNucSites` | bool | Stop simulation to view the nucleation sites that were set | true/false |


### LaserProperties
| Parameter Name | Value Type | Parameter Description | Unit |
| ------ | ------ | ------ | ------ | 
| `V_incident` | vector |Initial incidence vector for laser source | - |
| `HS_a` | scalar | Laser beam radius | m |
| `HS_bg` | scalar | Laser Offset in x axis | m |
| `HS_velocity` | scalar | Travel velocity of laser heat source in z direction | m/s |
| `HS_lg` | scalar | Laser Offset in z axis | m |
| `HS_Q` | scalar | Laser Power | W |
| `wavelength` | scalar | Wavelength of incident laser | m |
| `e_num_density` | scalar | electron number density of substrate |m^{-3} |
| `elec_resistivity` | scalar | Electrical resistivity of substrate | |

### g
| Parameter Name | Value Type | Parameter Description | Unit |
| ------ | ------ | ------ | ------ | 
| `dimensions` | OpenFOAM dimensions | Dimension set used to construct dimensionedVector | m/s |
| `value` | vector | Direction and value of gravitational acceleration | m/s2 |

### momentumTransport
| Parameter Name | Value Type | Parameter Description | Unit |
| ------ | ------ | ------ | ------ | 
| `simulationType` | Word | laminar?? | - |

### phaseProperties
| Parameter Name | Value Type | Parameter Description | Unit |
| ------ | ------ | ------ | ------ | 
| `phases` | List[Word] | Names of phases, a separate physicalProperties.[PHASE] dict will be read for each [PHASE]  | - |
| `sigma` | scalar | | |
| `dsigmadT` | scalar | | |
| `p0` | scalar | | |
| `Tvap` | scalar | | |
| `Mm` | scalar | | |
| `LatentHeatVap` | scalar | | |


### PhaseFieldProperties
| Parameter Name | Value Type | Parameter Description | Unit |
| ------ | ------ | ------ | ------ | 
| `grainNumThreshold` | scalar | | |
| `Aniso_strength` | scalar | | |
| `gamma_pf` | scalar | | |
| `sigma_P` | scalar | | |
| `sigma_G` | scalar | | |
| `lp` | scalar | | |
| `lg` | scalar | | |
| `deltaf_p` | scalar | | |
| `deltaf_p` | scalar | | |
| `D0` | scalar | | |
| `Qg` | scalar | | |
| `sigma_G` | scalar | | |
| `sigma_G` | scalar | | |
| `N_Seeds` | int | | |
| `Xmin` | scalar | | |
| `Xmax` | scalar | | |
| `Ymin` | scalar | | |
| `Ymax` | scalar | | |
| `Zmin` | scalar | | |
| `Zmax` | scalar | | |
| `sigma_G` | scalar | | |
| `sigma_G` | scalar | | |
| `write_ni_all` | bool | Set to true by default, to turn off the writing of INACTIVE n.i field, include `write_ni_all false;` | true/false |
| `write_ni_active` | bool | Set to true by default, to turn off the writing of ALL n.i fields, include `write_ni_active false;` as well as `write_ni_all false;`  | true/false |

### physicalProperties.[PHASE]
| Parameter Name | Value Type | Parameter Description | Unit |
| ------ | ------ | ------ | ------ | 
| `viscosityModel` | Word | | |
| `nu` | scalar | Dynamic viscosity | Pa.s |
| `rho` | scalar | Density | kg/m3 |
| `cp` | scalar | Specific heat capacity (constant pressure) | J/(kg.K) |
| `cpsolid` | scalar | Specific heat capacity...?? | J/(kg.K) |
| `kappa` | scalar | | |
| `kappasolid` | scalar | | |
| `Tsolidus` | scalar | Solidus temperature | K |
| `Tliquidus` | scalar | Liquidus temperature | K |
| `LatentHeat` | scalar | Latent heat of solidfification | J/kg |
| `beta` | scalar | Bouyancy?? | |


### Powder-Bed Fusion Example

.... text to go here ....

### Directional-Solidification Example

.... text to go here ....



## Algorithm

Initially the solver loads the mesh, reads in fields and boundary conditions, selects the turbulence model (if specified). The main solver loop is then initiated. First, the time step is
dynamically modified to ensure numerical stability. Next, the two-phase fluid mixture properties and turbulence quantities are updated. The discretized phase-fraction equation is then solved for a user-defined number of subtime steps (typically 3) using the multidimensional universal limiter with explicit solution solver [MULES](https://openfoam.org/release/2-3-0/multiphase/). This solver is included in the OpenFOAM library, and performs conservative solution of hyperbolic convective transport equations with defined bounds (0 and 1 for α1). Once the updated phase field is obtained, the program enters the pressure–velocity loop, in which p and u are corrected in an alternating fashion. In this loop T is also solved for, such that he buoyancy predictions are correct for the U and p fields. The process of correcting the pressure and velocity fields in sequence is known as pressure implicit with splitting of operators (PISO). In the OpenFOAM environment, PISO is repeated for multiple iterations at each time step. This process is referred to as merged PISO- semi-implicit method for pressure-linked equations (SIMPLE), or the pressure-velocity loop (PIMPLE) process, where SIMPLE is an iterative pressure–velocity solution algorithm. PIMPLE continues for a user specified number of iterations. 
The main solver loop iterates until program termination. A summary of the simulation algorithm is presented below:
* laserbeamFoam Simulation Algorithm Summary:
  * Initialize simulation data and mesh 
  * WHILE t<t_end DO
  * 1. Update delta_t for stability
  * 2. VOF equation sub-cycle
  * 3. Update interface location for heat source application
  * 4. Update fluid properties
  * 5. Ray-Tracing for Heat Source application at the surface
  * 6. PISO Loop
    * 1. Form u equation
    * 2. Energy Transport Loop
      * 1. Solve T equation
      * 2. Update fluid fraction field
      * 3. Re-evaluate source terms due to latent heat
    * 3. PISO
        * 1. Obtain and correct face fluxes
        * 2. Solve p-Poisson equation
        * 3. Correct u
  * 7. Nucleate New Order Parameters if undercooling is achieved
  * 8. Solve Phase-Field Equations
  * 9. Write Fields
  
There are no constraints on how the computational domain is discretised.

## License
OpenFoam, and by extension the laserbeamFoam application, is licensed free and open source only under the [GNU General Public Licence version 3](https://www.gnu.org/licenses/gpl-3.0.en.html). One reason for OpenFOAM’s popularity is that its users are granted the freedom to modify and redistribute the software and have a right of continued free use, within the terms of the GPL.

## Acknowledgements


## Citing This Work


## References

