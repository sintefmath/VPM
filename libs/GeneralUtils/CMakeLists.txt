project( GeneralUtils )
cmake_minimum_required( VERSION 3.0.2 )
find_path(MODHOT_ROOT_DIR /cmake/CMakeHeader.txt PATHS ../../)
include( "${MODHOT_ROOT_DIR}/cmake/CMakeHeader.txt")

add_library(GeneralUtils INTERFACE)
target_include_directories(GeneralUtils INTERFACE include/)
