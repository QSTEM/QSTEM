cmake_minimum_required(VERSION 2.8)

#project(qstem_experiments)

#find_package(Boost COMPONENTS system filesystem REQUIRED)

#include_directories(${Boost_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/libs "${FFTW3_INCLUDE_DIRS}")	

#FILE(GLOB LIBS_C_FILES "*.cpp")
#FILE(GLOB LIBS_H_FILES "*.hpp")

#set (qstem_expts_src ${LIBS_C_FILES} ${LIBS_H_FILES})

#add_library(qstem_experiments OBJECT ${qstem_expts_src})
#target_link_libraries(qstem_experiments ${Boost_LIBRARIES} qstem_potentials qstem_wavefunctions)

FILE(GLOB C_FILES "*.cpp")
FILE(GLOB H_FILES "*.hpp")

SOURCE_GROUP(experiments FILES ${C_FILES} ${H_FILES})


# Add files from this folder to the list of files for the shared library
set (LIBS_C_FILES ${LIBS_C_FILES} ${C_FILES} PARENT_SCOPE})
set (LIBS_C_FILES ${LIBS_H_FILES} ${H_FILES} PARENT_SCOPE})
