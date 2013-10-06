cmake_minimum_required(VERSION 2.8)

project(qstem_libs)

include_directories("${CMAKE_SOURCE_DIR}/libs" "${FFTW3_INCLUDE_DIRS}")	

FILE(GLOB STEM3_LIBS_C_FILES "${CMAKE_SOURCE_DIR}/libs/*.cpp")
FILE(GLOB STEM3_LIBS_H_FILES "${CMAKE_SOURCE_DIR}/libs/*.hpp")

set (qstem_libs_src ${STEM3_LIBS_C_FILES} ${STEM3_LIBS_H_FILES})
add_library(qstem_libs ${qstem_libs_src})

add_subdirectory(file_IO)
