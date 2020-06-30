project( NavierStokes_VPM_3D )
cmake_minimum_required( VERSION 2.8.11 )
find_path(MODHOT_ROOT_DIR /cmake/CMakeHeader.txt PATHS ../../)
include( "${MODHOT_ROOT_DIR}/cmake/CMakeHeader.txt")

add_modhot_lib(PoissonSolver3D)
add_modhot_lib(GeneralUtils)
add_modhot_lib(DiffusionEquationSolver3D)

file( GLOB sources "src/*.cpp" )
add_library( NavierStokes_VPM_3D ${sources} )
target_include_directories(NavierStokes_VPM_3D PUBLIC
    "./include/"
    )
target_link_libraries(NavierStokes_VPM_3D
    PoissonSolver3D
    GeneralUtils
    DiffusionEquationSolver3D
    )
