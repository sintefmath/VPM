add_library(GeneralUtils INTERFACE)
add_library(VPM::GeneralUtils ALIAS GeneralUtils)
target_include_directories(GeneralUtils INTERFACE include/)
