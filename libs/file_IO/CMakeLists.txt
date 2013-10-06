cmake_minimum_required(VERSION 2.8)

project(qstem_io_libs)

include_directories("${CMAKE_SOURCE_DIR}/libs" "${FFTW3_INCLUDE_DIRS}")	

FILE(GLOB QSTEM_LIB_HEADERS "${CMAKE_SOURCE_DIR}/libs/*.hpp")
FILE(GLOB STEM3_IO_LIBS_C_FILES "${CMAKE_SOURCE_DIR}/libs/file_IO/*.cpp")
FILE(GLOB STEM3_IO_LIBS_H_FILES "${CMAKE_SOURCE_DIR}/libs/file_IO/*.hpp")

set (qstem_io_libs_src ${STEM3_IO_LIBS_C_FILES} ${STEM3_IO_LIBS_H_FILES} ${QSTEM_LIB_HEADERS})
add_library(qstem_io_libs ${qstem_io_libs_src})
