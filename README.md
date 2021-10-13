# VPM
Vortex Particle Method for solving Navier Stokes in 2D and 3D

## Requirements

* [PETSc](https://petsc.org/release/)
* [pybind11](https://github.com/pybind/pybind11)
* A C++17 compatible compiler
* numpy

## Building

The project builds with cmake. It is probably a good idea to specify the following variables
* use cmake to set
    * ```-DCMAKE_CXX_COMPILER=/usr/bin/mpic++```
    * ```-DCMAKE_C_COMPILER=/usr/bin/mpicc```

Also make sure the install directory of pybind11 is in the ```CMAKE_PREFIX_PATH``` by specifying if necessary.

if you use something else than Ubuntu, you might need to:

* tweak ```FindPETSC.cmake``` in order for cmake to find PETSc
* tweak the openmpi path in the ```CMakeFile.txt``` of ```PoissonSolver2D```, etc.

## Running the Python scripts

We recommend that you create a virtualenvironment to install the dependencies (run the commands from the root of this repository):

    # Done once per setup
    python -m venv .venv
    # Done for every new terminal
    source .venv/bin/activate
    # Install packages (only needed once per setup)
    pip install -r requirements.txt

To run python scripts that depend on the ```pyVPM``` module, you need to make sure the folder ```build/pyVPM``` is in your ```PYTHONPATH```. From within the virtual environment activated above, you can run

    PYTHONPATH=$PYTHONPATH:/path/to/VPM/build/pyVPM python <name of script>

To run the ```plot_particles2d.py``` script, one could therefore run from ```apps/NavierStokes_2D_GeoRL/python_lib/```

    PYTHONPATH=$PYTHONPATH:/path/to/VPM/build/pyVPM python plot_particles2d.py




