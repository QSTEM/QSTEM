cmake_minimum_required(VERSION 2.8)

project(qstem_libs)

if(MINGW)
	include_directories("${CMAKE_SOURCE_DIR}/libs" "${CMAKE_SOURCE_DIR}/pthreads_w32" "${FFTW3_INCLUDE_DIRS}")
else(MINGW)
	include_directories("${CMAKE_SOURCE_DIR}/libs" "${FFTW3_INCLUDE_DIRS}")	
endif(MINGW)

set (fileio_LIB_SRC fileio_fftw3.h fileio_fftw3.cpp stemtypes_fftw3.h)
add_library(fileio STATIC ${fileio_LIB_SRC})

set (imagelib_LIB_SRC imagelib_fftw3.h imagelib_fftw3.cpp stemtypes_fftw3.h)
add_library(imagelib STATIC ${imagelib_LIB_SRC})

set (memory_LIB_SRC memory_fftw3.h memory_fftw3.cpp)
add_library(memory STATIC ${memory_LIB_SRC})

set (matrixlib_LIB_SRC matrixlib.h matrixlib.cpp stemtypes_fftw3.h)
add_library(matrixlib STATIC ${matrixlib_LIB_SRC})

set (readparams_LIB_SRC readparams.h readparams.cpp stemtypes_fftw3.h)
add_library(readparams STATIC ${readparams_LIB_SRC})

 