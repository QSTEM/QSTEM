cmake_minimum_required(VERSION 2.8)

project(qstem_libs)

find_package(Boost COMPONENTS system filesystem REQUIRED)

FIND_PACKAGE(HDF5 "1.8.5" REQUIRED)

# Make boost, fftw3, and our headers available to all projects

include_directories(${Boost_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/libs "${FFTW3_INCLUDE_DIRS}")
include_directories("${CMAKE_SOURCE_DIR}/libs" "${FFTW3_INCLUDE_DIRS}" "${HDF5_INCLUDE_DIR}" )	

FILE(GLOB LIBS_C_FILES "*.cpp")
FILE(GLOB LIBS_H_FILES "*.hpp")
# Include all of the IO functionality with the main library
FILE(GLOB CONFIG_C_FILES "config_IO/*.cpp")
FILE(GLOB CONFIG_H_FILES "config_IO/*.hpp")
FILE(GLOB DATA_C_FILES "data_IO/*.cpp")
FILE(GLOB DATA_H_FILES "data_IO/*.hpp")
FILE(GLOB STRUCTURE_C_FILES "structure_IO/*.cpp")
FILE(GLOB STRUCTURE_H_FILES "structure_IO/*.hpp")
FILE(GLOB POTENTIAL_C_FILES "potentials/*.cpp")
FILE(GLOB POTENTIAL_H_FILES "potentials/*.hpp")
FILE(GLOB EXPERIMENT_C_FILES "experiments/*.cpp")
FILE(GLOB EXPERIMENT_H_FILES "experiments/*.hpp")

set (qstem_libs_src ${LIBS_C_FILES} ${LIBS_H_FILES})
set (qstem_libs_src ${qstem_libs_src} ${CONFIG_C_FILES} ${CONFIG_H_FILES} )
set (qstem_libs_src ${qstem_libs_src} ${DATA_C_FILES} ${DATA_H_FILES} )
set (qstem_libs_src ${qstem_libs_src} ${STRUCTURE_C_FILES} ${STRUCTURE_H_FILES} )
set (qstem_libs_src ${qstem_libs_src} ${POTENTIAL_C_FILES} ${POTENTIAL_H_FILES} )
set (qstem_libs_src ${qstem_libs_src} ${EXPERIMENT_C_FILES} ${EXPERIMENT_H_FILES} )

add_library(qstem_libs ${qstem_libs_src})
target_link_libraries(qstem_libs ${Boost_LIBRARIES})

# Also build IO functionality separately for external use (??)
#add_subdirectory(config_IO)
#add_subdirectory(data_IO)
#add_subdirectory(structure_IO)

