project( PoissonSolver3D )
cmake_minimum_required( VERSION 2.8.11 )
find_path(MODHOT_ROOT_DIR /cmake/CMakeHeader.txt PATHS ../../)
include( "${MODHOT_ROOT_DIR}/cmake/CMakeHeader.txt")

FIND_PACKAGE( PETSc REQUIRED ) 

file( GLOB sources "src/*.cpp" )
add_library( PoissonSolver3D ${sources} )
target_include_directories(PoissonSolver3D PUBLIC
    "./include/"
    "${PETSC_DIR}/include"
    "${PETSC_CMAKE_ARCH}/include"
    )
target_link_libraries(PoissonSolver3D
    ${PETSC_LIBRARIES}
    )
