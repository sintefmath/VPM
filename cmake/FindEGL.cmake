# Find EGL
#
# EGL_INCLUDE_DIR
# EGL_LIBRARY
# EGL_FOUND

find_path(EGL_INCLUDE_DIR NAMES EGL/egl.h)

set(EGL_NAMES ${EGL_NAMES} egl EGL)
find_library(EGL_LIBRARY NAMES ${EGL_NAMES})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EGL DEFAULT_MSG EGL_LIBRARY EGL_INCLUDE_DIR)

add_definitions( -D_MODHOT_EGL_ )

set(EGL_INCLUDE_DIRS ${EGL_INCLUDE_DIR})
set(EGL_LIBRARIES ${EGL_LIBRARY})

mark_as_advanced(EGL_INCLUDE_DIR EGL_LIBRARY)
