project( NavierStokes_VPM_2D )
cmake_minimum_required( VERSION 2.8.11 )
find_path(MODHOT_ROOT_DIR /cmake/CMakeHeader.txt PATHS ../../)
include( "${MODHOT_ROOT_DIR}/cmake/CMakeHeader.txt")

add_modhot_lib(PoissonSolver2D)
add_modhot_lib(GeneralUtils)
add_modhot_lib(DiffusionEquationSolver2D)

file( GLOB sources "src/*.cpp" )
add_library( NavierStokes_VPM_2D ${sources} )
target_include_directories(NavierStokes_VPM_2D PUBLIC
    "./include/"
    )
target_link_libraries(NavierStokes_VPM_2D
    PoissonSolver2D
    GeneralUtils
    DiffusionEquationSolver2D
    )
