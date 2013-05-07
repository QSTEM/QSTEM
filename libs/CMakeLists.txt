cmake_minimum_required(VERSION 2.8)

project(qstem_libs)

if(MINGW)
	include_directories("${CMAKE_SOURCE_DIR}/libs" "${CMAKE_SOURCE_DIR}/pthreads_w32" "${FFTW3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}")
else(MINGW)
	include_directories("${CMAKE_SOURCE_DIR}/libs" "${FFTW3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}")	
endif(MINGW)

FILE(GLOB STEM3_LIBS_C_FILES "${CMAKE_SOURCE_DIR}/libs/*.cpp")
FILE(GLOB STEM3_LIBS_H_FILES "${CMAKE_SOURCE_DIR}/libs/*.h")

set (qstem_libs_src ${STEM3_LIBS_C_FILES} ${STEM3_LIBS_H_FILES})
add_library(qstem_libs STATIC ${qstem_libs_src})

 