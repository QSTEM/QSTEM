cmake_minimum_required(VERSION 2.8)

project(stem3)		

include_directories(${CMAKE_SOURCE_DIR}/libs)

FILE(GLOB QSTEM_LIB_HEADERS "${CMAKE_SOURCE_DIR}/libs/*.hpp")

FILE(GLOB STEM3_C_FILES "${CMAKE_SOURCE_DIR}/stem3/*.cpp")
FILE(GLOB STEM3_H_FILES "${CMAKE_SOURCE_DIR}/stem3/*.hpp")

add_executable(stem3 ${STEM3_C_FILES} ${STEM3_H_FILES} ${QSTEM_LIB_HEADERS})
# m is libm - math libraries on Unix systems
target_link_libraries(stem3 qstem_libs ${FFTW3_LIBS} ${FFTW3F_LIBS} ${M_LIB} ${Boost_LIBRARIES})
 
if(OPENMP)
	SET_TARGET_PROPERTIES(stem3 PROPERTIES COMPILE_FLAGS "${OpenMP_C_FLAGS}" LINK_FLAGS  "${OpenMP_C_FLAGS}")
endif(OPENMP)
