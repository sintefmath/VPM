FIND_PACKAGE( PETSc REQUIRED ) 

file( GLOB sources "src/*.cpp" )
add_library(PoissonSolver3D ${sources} )
target_include_directories(PoissonSolver3D PUBLIC
    "./include/"
    "${PETSC_DIR}/include"
    "${PETSC_CMAKE_ARCH}/include"
    )
target_link_libraries(PoissonSolver3D PUBLIC
    ${PETSC_LIBRARIES}
    )

add_library(VPM::PoissonSolver3D ALIAS PoissonSolver3D)
