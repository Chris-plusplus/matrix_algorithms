include_guard()

include("${PROJECT_SOURCE_DIR}/cmake/conan.cmake")

add_library(${PROJECT_NAME})

# find source files
file(GLOB_RECURSE MATRIX_SOURCE src/**.cpp)
target_sources(${PROJECT_NAME} PRIVATE ${MATRIX_SOURCE})
target_include_directories(${PROJECT_NAME} PUBLIC include)

# link conan libraries
target_link_libraries(${PROJECT_NAME} PUBLIC ${MATRIX_LIBRARIES})
