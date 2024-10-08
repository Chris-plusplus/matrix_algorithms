include_guard()

# check build type
if(
    "${CMAKE_CONFIGURATION_TYPES}" STREQUAL Debug
    OR "${CMAKE_BUILD_TYPE}" STREQUAL Debug
    OR "${CMAKE_CONFIGURATION_TYPES}" STREQUAL "RelWithDebInfo"
    OR "${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo"
)
    set(MATRIX_DEBUG TRUE)
    set(MATRIX_RELEASE FALSE)

    set(MATRIX_BUILD_TYPE "Debug")
    
    add_compile_definitions(MATRIX_DEBUG=1)
    add_compile_definitions(MATRIX_RELEASE=0)
else()
    set(MATRIX_DEBUG FALSE)
    set(MATRIX_RELEASE TRUE)
    
    set(MATRIX_BUILD_TYPE "Release")
    
    add_compile_definitions(MATRIX_DEBUG=0)
    add_compile_definitions(MATRIX_RELEASE=1)
endif()

message(STATUS "Build mode: ${MATRIX_BUILD_TYPE}")

set(CMAKE_BUILD_TYPE ${MATRIX_BUILD_TYPE})
