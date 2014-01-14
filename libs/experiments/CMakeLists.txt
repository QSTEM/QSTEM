cmake_minimum_required(VERSION 2.8)

project(qstem_experiments)

find_package(Boost COMPONENTS system filesystem REQUIRED)

# Make boost, fftw3, and our headers available to all projects

include_directories(${Boost_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/libs "${FFTW3_INCLUDE_DIRS}")
include_directories("${CMAKE_SOURCE_DIR}/libs" "${FFTW3_INCLUDE_DIRS}" )	

FILE(GLOB LIBS_C_FILES "*.cpp")
FILE(GLOB LIBS_H_FILES "*.hpp")

set (qstem_expts_src ${LIBS_C_FILES} ${LIBS_H_FILES})

add_library(qstem_experiments ${qstem_expts_src})
target_link_libraries(qstem_experiments ${Boost_LIBRARIES})
target_link_libraries(qstem_experiments qstem_potentials)
target_link_libraries(qstem_experiments qstem_wavefunctions)


# Also build IO functionality separately for external use (??)
#add_subdirectory(config_IO)
#add_subdirectory(data_IO)
#add_subdirectory(structure_IO)

