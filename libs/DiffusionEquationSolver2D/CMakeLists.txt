FIND_PACKAGE( PETSc REQUIRED ) 

file( GLOB sources "src/*.cpp" )
add_library( DiffusionEquationSolver2D ${sources} )
target_include_directories( DiffusionEquationSolver2D PUBLIC
    "./include/"
    "${PETSC_DIR}/include"
    "${PETSC_CMAKE_ARCH}/include"
    )
target_link_libraries( DiffusionEquationSolver2D
    ${PETSC_LIBRARIES}
    )

add_library(VPM::DiffusionEquationSolver2D ALIAS DiffusionEquationSolver2D)
