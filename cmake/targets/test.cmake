include_guard()

#include("${PROJECT_SOURCE_DIR}/cmake/library.cmake")

add_executable(test_bin programs/test.cpp)
#target_link_libraries(test_bin PUBLIC ${PROJECT_NAME})
