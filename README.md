# VPM
Vortex Particle Method for solving Navier Stokes in 2D and 3D

* install petsc

* use cmake to set
    * CMAKE_CXX_COMPILER              */usr/bin/mpic++
    * CMAKE_C_COMPILER                */usr/bin/mpicc

if you use something else than Ubuntu, you might need to:
    * you might need to tweak FindPETSC.cmake in order for cmake to find petsc
    * as well as the openmpi path in the CMakeFile.txt of PoissonSolver2D

