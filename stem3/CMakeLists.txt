cmake_minimum_required(VERSION 2.8)

project(stem3)		

if(MINGW)
	include_directories("${CMAKE_SOURCE_DIR}/libs" "${CMAKE_SOURCE_DIR}/pthreads_w32" "${FFTW3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}")
	find_library(PTHREAD_LIB libpthread "${CMAKE_SOURCE_DIR}/pthreads_w32")
	add_definitions(-DPTW32_STATIC_LIB)
else(MINGW)
	include_directories("${CMAKE_SOURCE_DIR}/libs" "${FFTW3_INCLUDE_DIRS}" "${Boost_INCLUDE_DIRS}")
endif(MINGW)

FILE(GLOB QSTEM_LIB_HEADERS "${CMAKE_SOURCE_DIR}/libs/*.h")

FILE(GLOB STEM3_C_FILES "${CMAKE_SOURCE_DIR}/stem3/*.cpp")
FILE(GLOB STEM3_H_FILES "${CMAKE_SOURCE_DIR}/stem3/*.h")

# m is libm - math libraries on Unix systems

add_executable(stem3 ${STEM3_C_FILES} ${STEM3_H_FILES} ${QSTEM_LIB_HEADERS})
if(MINGW)
	target_link_libraries(stem3 qstem_libs ${FFTW3_LIBS} ${FFTW3F_LIBS}
                                ${PTHREAD_LIB} ${M_LIB})
else(MINGW)
	target_link_libraries(stem3 qstem_libs	${FFTW3_LIBS} ${FFTW3F_LIBS} ${M_LIB})
endif(MINGW)
 
 if(OPENMP)
	SET_TARGET_PROPERTIES(stem3 PROPERTIES COMPILE_FLAGS "${OpenMP_C_FLAGS}"
										  LINK_FLAGS "${OpenMP_C_FLAGS}")

endif(OPENMP)
