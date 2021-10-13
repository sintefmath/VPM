add_executable(navier_stokes_3d "./src/App.cpp" )
target_link_libraries( navier_stokes_3d
    VPM::NavierStokes_VPM_3D
    )
