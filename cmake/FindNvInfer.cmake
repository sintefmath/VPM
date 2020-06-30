# Find NVInfer

find_library(NVINFER_LIBRARY NAMES libnvinfer.so Engine.lib nvinfer.lib
    HINTS
    ${NVINFER_ROOT}/lib/x64
    )
find_library(NVONNXPARSER_LIBRARY NAMES nvonnxparser.lib
    HINTS
    ${NVINFER_ROOT}/lib/x64
    )
find_library(NVINFER_PLUGIN_LIBRARY NAMES nvinfer_plugin.lib
    HINTS
    ${NVINFER_ROOT}/lib/x64
    )

set(NVINFER_LIBRARIES ${NVINFER_LIBRARY} ${NVONNXPARSER_LIBRARY} ${NVINFER_PLUGIN_LIBRARY})
set(NVINFER_INCLUDE_DIRS ${NVINFER_INCLUDE_DIR})

mark_as_advanced(NVINFER_LIBRARIES NVINFER_INCLUDE_DIRS)
