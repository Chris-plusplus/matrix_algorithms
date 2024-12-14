cmake_minimum_required(VERSION 3.19)

project(
        matrix_algorithms
        VERSION 0.1.0
        LANGUAGES CXX C
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")


include("cmake/build_mode.cmake")
include("cmake/os.cmake")
include("cmake/msvc.cmake")
include("cmake/defines.cmake")

# uncomment below to force conan reinstallation
set(MATRIX_FORCE_CONAN_INSTALL TRUE)
include("cmake/conan.cmake")

include("cmake/library.cmake")

# CMake targets
file(GLOB_RECURSE MATRIX_TARGETS cmake/targets/*.cmake)
foreach (MATRIX_TARGET ${MATRIX_TARGETS})
    include(${MATRIX_TARGET})
endforeach ()
