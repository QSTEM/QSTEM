cmake_minimum_required(VERSION 2.8)

project(qstem_structure_lib)

include_directories("${CMAKE_SOURCE_DIR}/libs" "${FFTW3_INCLUDE_DIRS}")	

FILE(GLOB QSTEM_LIB_HEADERS "${CMAKE_SOURCE_DIR}/libs/*.hpp")
FILE(GLOB STEM3_STRUCTURE_LIBS_C_FILES "*.cpp")
FILE(GLOB STEM3_STRUCTURE_LIBS_H_FILES "*.hpp")

set (qstem_structure_libs_src ${STEM3_STRUCTURE_LIBS_C_FILES} ${STEM3_STRUCTURE_LIBS_H_FILES} ${QSTEM_LIB_HEADERS})
add_library(qstem_structure_lib ${qstem_structure_libs_src})
