project( PoissonSolver2D )
cmake_minimum_required( VERSION 2.8.11 )
find_path(MODHOT_ROOT_DIR /cmake/CMakeHeader.txt PATHS ../../)
include( "${MODHOT_ROOT_DIR}/cmake/CMakeHeader.txt")

FIND_PACKAGE( PETSc REQUIRED ) 

file( GLOB sources "src/*.cpp" )
add_library( PoissonSolver2D ${sources} )
target_include_directories(PoissonSolver2D PUBLIC
    "./include/"
    "${PETSC_DIR}/include"
    "${PETSC_INCLUDE_DIR}/include"
    "${PETSC_CMAKE_ARCH}/include"
    "/usr/lib/x86_64-linux-gnu/openmpi/include/"
    )
target_link_libraries(PoissonSolver2D
    ${PETSC_LIBRARIES}
    )
