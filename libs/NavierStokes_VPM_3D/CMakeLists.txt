file( GLOB sources "src/*.cpp" )

add_library( NavierStokes_VPM_3D ${sources} )
target_include_directories(NavierStokes_VPM_3D PUBLIC
    "./include/"
    )
target_link_libraries(NavierStokes_VPM_3D PUBLIC
    VPM::PoissonSolver3D
    VPM::GeneralUtils
    VPM::DiffusionEquationSolver3D
    )

add_library(VPM::NavierStokes_VPM_3D ALIAS NavierStokes_VPM_3D)
