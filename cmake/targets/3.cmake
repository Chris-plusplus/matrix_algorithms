include_guard()

include("${PROJECT_SOURCE_DIR}/cmake/library.cmake")

add_executable(3_bin programs/3.cpp)
target_link_libraries(3_bin PUBLIC ${PROJECT_NAME})
