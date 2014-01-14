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
FILE(GLOB CONFIG_H_FILES "config_IO/*.hpp")
FILE(GLOB DATA_H_FILES "data_IO/*.hpp")
FILE(GLOB STRUCTURE_H_FILES "structure_IO/*.hpp")
FILE(GLOB POTENTIAL_H_FILES "potentials/*.hpp")
FILE(GLOB EXPERIMENT_H_FILES "experiments/*.hpp")

set (qstem_libs_src ${LIBS_C_FILES} ${LIBS_H_FILES})
set (qstem_libs_src ${qstem_libs_src} ${CONFIG_H_FILES} )
set (qstem_libs_src ${qstem_libs_src} ${DATA_H_FILES} )
set (qstem_libs_src ${qstem_libs_src} ${STRUCTURE_H_FILES} )
set (qstem_libs_src ${qstem_libs_src} ${POTENTIAL_H_FILES} )
set (qstem_libs_src ${qstem_libs_src} ${EXPERIMENT_H_FILES} )

# defines qstem_config_lib (config readers)
add_subdirectory(config_IO)

# defines qstem_data_lib (I/O for image data)
add_subdirectory(data_IO)

#defines qstem_structure_lib
add_subdirectory(structure_IO)

# defines qstem_wavefunctions_lib
add_subdirectory(wavefunctions)

# defines qstem_potentials_lib
add_subdirectory(potentials)

# defines qstem_experiments_lib
add_subdirectory(experiments)


add_library(qstem ${qstem_libs_src})
target_link_libraries(qstem ${Boost_LIBRARIES})
target_link_libraries(qstem qstem_config)
target_link_libraries(qstem qstem_data)
target_link_libraries(qstem qstem_structure)
target_link_libraries(qstem qstem_wavefunctions)
target_link_libraries(qstem qstem_potentials)
target_link_libraries(qstem qstem_experiments)
