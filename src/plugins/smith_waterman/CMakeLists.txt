set(UGENE_PLUGIN_ID smith_waterman)
set(UGENE_PLUGIN_NAME "Smith-Waterman")

file(GLOB_RECURSE SRCS src/*.cpp src/*.c src/*.h)

# SSE2 sources must not be compiled for Apple Silicon (or other non-x86).
# Use CMAKE_OSX_ARCHITECTURES on macOS: CMAKE_SYSTEM_PROCESSOR stays arm64 even
# when an Intel-only build is requested.
if (APPLE)
    if (NOT CMAKE_OSX_ARCHITECTURES MATCHES "x86_64")
        list(REMOVE_ITEM SRCS "${CMAKE_CURRENT_SOURCE_DIR}/src/SmithWatermanAlgorithmSSE2.h")
        list(REMOVE_ITEM SRCS "${CMAKE_CURRENT_SOURCE_DIR}/src/SmithWatermanAlgorithmSSE2.cpp")
    endif ()
elseif (NOT CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
    list(REMOVE_ITEM SRCS "${CMAKE_CURRENT_SOURCE_DIR}/src/SmithWatermanAlgorithmSSE2.h")
    list(REMOVE_ITEM SRCS "${CMAKE_CURRENT_SOURCE_DIR}/src/SmithWatermanAlgorithmSSE2.cpp")
endif ()

include(${CMAKE_SOURCE_DIR}/cmake/Plugin.cmake)
