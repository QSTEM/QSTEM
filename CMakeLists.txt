cmake_minimum_required(VERSION 2.8)

set(CMAKE_USER_MAKE_RULES_OVERRIDE
	${CMAKE_SOURCE_DIR}/cmake/Modules/c_flag_overrides.cmake)
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX
	${CMAKE_SOURCE_DIR}/cmake/Modules/cxx_flag_overrides.cmake)

project(qstem)
set (stem3_VERSION_MAJOR 2)
set (stem3_VERSION_MINOR 23)

set (EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin/")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
					
find_package(fftw3 REQUIRED)

SET(BOOST_MIN_VERSION "1.38.0")
FIND_PACKAGE(Boost ${BOOST_MIN_VERSION} REQUIRED)

IF (Boost_FOUND)
	message(STATUS "Found Boost.  Ignore message that it is not found.")
    INCLUDE_DIRECTORIES(${BOOST_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()

set (M_LIB "")

if (UNIX)
	set (M_LIB "m")
endif(UNIX)


if (NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Debug")
	set(CMAKE_BUILD_TYPE "Debug")
endif(NOT CMAKE_BUILD_TYPE)

OPTION( OPENMP "Set to ON to enable parallel execution using OpenMP" ON )

if(OPENMP)
	find_package(OpenMP REQUIRED)
endif(OPENMP)

#set(CMAKE_C_FLAGS_DEBUG "-ggdb -static")
#set(CMAKE_CXX_FLAGS_DEBUG "-ggdb -static")
#set(CMAKE_C_FLAGS_RELEASE "-O3 -static -Wall -march=native")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -static -Wall -march=native")

# adds the libraries
add_subdirectory(libs)
add_subdirectory(stem3)
add_subdirectory(gbmaker)
add_subdirectory(qscRg12)
 
