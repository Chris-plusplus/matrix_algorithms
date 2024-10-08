include_guard()

# check for OS
if(WIN32)
    set(MATRIX_WINDOWS TRUE)
    set(MATRIX_LINUX FALSE)
elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    set(MATRIX_WINDOWS FALSE)
    set(MATRIX_LINUX TRUE)
endif()

# set user home directory
if(MATRIX_WINDOWS)    
    add_compile_definitions(MATRIX_WINDOWS=1)
    add_compile_definitions(MATRIX_LINUX=0)
    add_compile_definitions(MATRIX_UNIX=0)

    message(STATUS "OS: Windows")
elseif(MATRIX_LINUX)    
    add_compile_definitions(MATRIX_WINDOWS=0)
    add_compile_definitions(MATRIX_LINUX=1)
    add_compile_definitions(MATRIX_UNIX=1)

    message(STATUS "OS: Linux")
endif()
