file( GLOB sources "src/*.cpp" )
add_library(NavierStokes_VPM_2D SHARED ${sources} )
target_include_directories(NavierStokes_VPM_2D PUBLIC
    "./include/"
    )
target_link_libraries(NavierStokes_VPM_2D PUBLIC
    VPM::PoissonSolver2D
    VPM::GeneralUtils
    VPM::DiffusionEquationSolver2D
    VPM::BBFMM2D
    )
target_compile_definitions(NavierStokes_VPM_2D PRIVATE -D_USE_BBFMM_POSTYPE_)
set_target_properties(NavierStokes_VPM_2D PROPERTIES CXX_STANDARD 17 )
add_library(VPM::NavierStokes_VPM_2D ALIAS NavierStokes_VPM_2D)
