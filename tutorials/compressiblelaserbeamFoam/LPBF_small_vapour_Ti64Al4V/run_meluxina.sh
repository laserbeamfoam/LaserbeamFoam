#!/bin/bash -l
## This file is called `LPBF_tooth_1.sh`
#SBATCH --time=40:00:00
#SBATCH --account=p201141
#SBATCH --partition=cpu
#SBATCH --qos=default
#SBATCH -N 1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=cpetar112@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

cd $SLURM_SUBMIT_DIR

module load env/release/2022.1
module load env/staging/2022.1
#module load Arm-Forge/22.0.4-GCC-11.3.0
#source /project/home/p201141/Petar/OpenFOAM-v2506/etc/bashrc
source /project/home/p201141/Software/OpenFOAM-v2512/etc/bashrc
module load OpenMPI/4.1.4-GCC-11.3.0

#mpirun -n 65 laserMeltFoam_new_ray_tracing -parallel &> log.laserMeltFoam
mpirun -n 16 compressibleLaserbeamFoamUPDATED -parallel &> log.compressibleLaserbeamFoamUPDATED

#setSolidFraction -subDivisions 10

#laserMeltFoam_adaptive_dt_v4
#laserMeltFoam

#zip -r LPBF_small_vapour.zip LPBF_small_vapour > log
