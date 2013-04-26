cmake_minimum_required(VERSION 2.8)

project(gbmaker)

if(MINGW)
	include_directories("${CMAKE_SOURCE_DIR}/libs" "${CMAKE_SOURCE_DIR}/pthreads_w32" "${FFTW3_INCLUDE_DIRS}")
else(MINGW)
	include_directories("${CMAKE_SOURCE_DIR}/libs" "${FFTW3_INCLUDE_DIRS}")	
endif(MINGW)

add_executable(gbmaker gbmaker.cpp)
target_link_libraries(gbmaker qstem_libs ${FFTW3_LIBS}	${FFTW3F_LIBS} ${M_LIB})
 
