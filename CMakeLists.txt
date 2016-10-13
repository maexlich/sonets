cmake_minimum_required(VERSION 3.6)
project(sonets_cmake)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(GSL REQUIRED)

if(NOT GSL_FOUND)
    message(ERROR Sonets requires GSL to compile)
endif()

# indicate we wish for the binaries to be written to the bin directory
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

# build items in subdirectories
add_subdirectory(src)
