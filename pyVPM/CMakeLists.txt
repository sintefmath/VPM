find_package(Python COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG)

pybind11_add_module(pyVPM MODULE src/PyVPM.cpp)
target_link_libraries(pyVPM PUBLIC Python::Module VPM::NavierStokes_VPM_2D)
set_target_properties(pyVPM PROPERTIES CXX_STANDARD 17 )
