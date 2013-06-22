cmake_minimum_required(VERSION 2.8)

set(CMAKE_USER_MAKE_RULES_OVERRIDE
	${CMAKE_SOURCE_DIR}/cmake/Modules/c_flag_overrides.cmake)
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX
	${CMAKE_SOURCE_DIR}/cmake/Modules/cxx_flag_overrides.cmake)

project(qstem)
set (stem3_VERSION_MAJOR 2)
set (stem3_VERSION_MINOR 30)

#set (EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin/")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
					
find_package(fftw3 REQUIRED)

find_package(Boost 1.35.0)

# Make boost, fftw3, and our headers available to all projects
include_directories(${Boost_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/libs "${FFTW3_INCLUDE_DIRS}")

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

# Squelch Visual studio's warnings about insecure functions - will replace these over time, but must maintain Linux compatibility.
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

# adds the libraries
add_subdirectory(libs)
add_subdirectory(stem3)
add_subdirectory(gbmaker)
add_subdirectory(qscRg12)
add_subdirectory(tests)

enable_testing ()
add_test (NAME TestLibs COMMAND test_libs)
#add_test (NAME TestStem3 COMMAND test_stem3)
#add_test (NAME TestGBMaker COMMAND test_gbmaker)