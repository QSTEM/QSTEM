cmake_minimum_required(VERSION 2.8)

project(qstem_libs)

if(CMAKE_COMPILER_IS_GNUCXX)
    include(CodeCoverage)
    setup_target_for_coverage(${PROJECT_NAME}_coverage ${PROJECT_TEST_NAME} coverage)
endif()

include_directories("${CMAKE_SOURCE_DIR}/libs" "${FFTW3_INCLUDE_DIRS}")	

FILE(GLOB STEM3_LIBS_C_FILES "${CMAKE_SOURCE_DIR}/libs/*.cpp")
FILE(GLOB STEM3_LIBS_H_FILES "${CMAKE_SOURCE_DIR}/libs/*.h")

set (qstem_libs_src ${STEM3_LIBS_C_FILES} ${STEM3_LIBS_H_FILES})
add_library(qstem_libs ${qstem_libs_src})
