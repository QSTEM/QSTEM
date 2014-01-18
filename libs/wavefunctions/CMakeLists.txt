cmake_minimum_required(VERSION 2.8)

#project(qstem_wavefunctions)

#include_directories("${CMAKE_SOURCE_DIR}/libs" "${FFTW3_INCLUDE_DIRS}")	

#FILE(GLOB QSTEM_LIB_HEADERS "${CMAKE_SOURCE_DIR}/libs/*.hpp")
#FILE(GLOB STEM3_WAVEFUNCTIONS_LIBS_C_FILES "*.cpp")
#FILE(GLOB STEM3_WAVEFUNCTIONS_LIBS_H_FILES "*.hpp")

#set (qstem_wavefunctions_libs_src ${STEM3_WAVEFUNCTIONS_LIBS_C_FILES} ${STEM3_WAVEFUNCTIONS_LIBS_H_FILES} ${QSTEM_LIB_HEADERS})
#add_library(qstem_wavefunctions OBJECT ${qstem_wavefunctions_libs_src})

FILE(GLOB C_FILES "*.cpp")
FILE(GLOB H_FILES "*.hpp")

SOURCE_GROUP(wavefunctions FILES ${C_FILES} ${H_FILES})

# Add files from this folder to the list of files for the shared library
set (LIBS_C_FILES ${LIBS_C_FILES} ${C_FILES} PARENT_SCOPE})
set (LIBS_C_FILES ${LIBS_H_FILES} ${H_FILES} PARENT_SCOPE})
