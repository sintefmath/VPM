find_package(MPI REQUIRED)
add_executable(navier_stokes_2d_georl "./src/App.cpp" )
target_link_libraries( navier_stokes_2d_georl
    VPM::NavierStokes_VPM_2D MPI::MPI_C
    )

add_executable(init "./src/Init.cpp" )
target_link_libraries( init
    VPM::NavierStokes_VPM_2D MPI::MPI_C
    )
