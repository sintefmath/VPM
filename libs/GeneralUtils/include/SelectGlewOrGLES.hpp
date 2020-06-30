#pragma once

// Include OpenGL ES headers instead of GLEW if running under Emscripten or Android.

#ifdef MODHOT_GL
    #include <GL/glew.h>
#else
    #include <GLES3/gl3.h>
    #define GL_GLEXT_PROTOTYPES
    #include <GLES2/gl2ext.h>
#endif
