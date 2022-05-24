find_package(MPI REQUIRED)

add_executable(navier_stokes_2d "./src/App.cpp" )
target_link_libraries( navier_stokes_2d
    VPM::NavierStokes_VPM_2D MPI::MPI_C
    )
