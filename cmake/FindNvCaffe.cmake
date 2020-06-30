# Find NvCaffe

find_library(NVCAFFE_LIBRARY NAMES libnvcaffe_parser.so CaffeParser.lib nvparsers.lib
    HINTS
    ${NVCAFFEE_ROOT}/lib/x64
    )
find_path(NVCAFFE_INCLUDE_DIR NvCaffeParser.h
    HINTS
    ${NVCAFFEE_ROOT}/include
    )

set(NVCAFFE_LIBRARIES ${NVCAFFE_LIBRARY})
set(NVCAFFE_INCLUDE_DIRS ${NVCAFFE_INCLUDE_DIR})

mark_as_advanced(NVCAFFE_LIBRARIES NVCAFFE_INCLUDE_DIRS)
