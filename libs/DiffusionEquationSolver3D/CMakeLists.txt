FIND_PACKAGE( PETSc REQUIRED ) 

file( GLOB sources "src/*.cpp" )
add_library(DiffusionEquationSolver3D ${sources} )
target_include_directories(DiffusionEquationSolver3D PUBLIC
    "./include/"
    "${PETSC_DIR}/include"
    "${PETSC_CMAKE_ARCH}/include"
    )
target_link_libraries(DiffusionEquationSolver3D PUBLIC
    ${PETSC_LIBRARIES}
    )

add_library(VPM::DiffusionEquationSolver3D ALIAS DiffusionEquationSolver3D)
