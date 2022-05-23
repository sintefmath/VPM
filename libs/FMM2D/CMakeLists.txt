find_package (Eigen3 3.3 REQUIRED NO_MODULE)
file( GLOB sources "BBFMM2D/src/*.cpp" )
add_library( BBFMM2D ${sources} )
target_include_directories( BBFMM2D PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/BBFMM2D/header/"
    )
target_link_libraries( BBFMM2D Eigen3::Eigen)

add_library(VPM::BBFMM2D ALIAS BBFMM2D)
