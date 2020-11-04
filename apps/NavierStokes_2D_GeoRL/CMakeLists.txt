project( NavierStokes_2D )
cmake_minimum_required( VERSION 2.8.11 )
find_path(MODHOT_ROOT_DIR /cmake/CMakeHeader.txt PATHS ../../)
include( "${MODHOT_ROOT_DIR}/cmake/CMakeHeader.txt")

find_package(MPI REQUIRED)
add_modhot_lib(NavierStokes_VPM_2D)
add_executable(app "./src/App.cpp" )
target_link_libraries( app
    NavierStokes_VPM_2D
    )

add_executable(init "./src/Init.cpp" )
target_link_libraries( init
    NavierStokes_VPM_2D MPI::MPI_C
    )
