cmake_minimum_required(VERSION 2.8)

#include_directories("${CMAKE_SOURCE_DIR}/libs" "${FFTW3_INCLUDE_DIRS}")	

#FILE(GLOB QSTEM_LIB_HEADERS "${CMAKE_SOURCE_DIR}/libs/*.hpp")
FILE(GLOB C_FILES "*.cpp")
FILE(GLOB H_FILES "*.hpp")

SOURCE_GROUP(potentials FILES ${C_FILES} ${H_FILES})

# Add files from this folder to the list of files for the shared library
set (LIBS_C_FILES ${LIBS_C_FILES} ${C_FILES} PARENT_SCOPE})

#set (qstem_potential_libs_src ${C_FILES} ${H_FILES} ${QSTEM_LIB_HEADERS})
#add_library(qstem_potentials OBJECT ${qstem_potential_libs_src})
#add_library(qstem_potentials ${qstem_potential_libs_src})

#set(QSTEM_LIBS ${QSTEM_LIBS} qstem_potentials PARENT_SCOPE)
