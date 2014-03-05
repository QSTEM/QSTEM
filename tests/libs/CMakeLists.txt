cmake_minimum_required(VERSION 2.8)

project(tests)

add_definitions (-DBOOST_TEST_DYN_LINK)
add_definitions(-DBOOST_ALL_NO_LIB)  # tell the compiler to undefine this boost macro

set(Boost_USE_STATIC_LIBS        OFF)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     OFF)
set(BOOST_ALL_DYN_LINK           ON)   # force dynamic linking for all libraries

find_package(Boost 1.35.0 COMPONENTS unit_test_framework REQUIRED)

include_directories("${CMAKE_SOURCE_DIR}/libs" "${FFTW3_INCLUDE_DIRS}" ${Boost_INCLUDE_DIRS})

set(DIRS config_IO wavefunctions data_IO potentials experiments structure_IO)

FOREACH(dirname ${DIRS})
	add_subdirectory(${dirname})
ENDFOREACH(dirname ${DIRS})

FILE(GLOB TEST_FILES "*.cpp")
FOREACH(TEST_FILE ${TEST_FILES})
  GET_FILENAME_COMPONENT(filename ${TEST_FILE} NAME_WE)
  add_executable(${filename} ${TEST_FILE})
  target_link_libraries(${filename} qstem ${FFTW3_LIBS} ${FFTW3F_LIBS} ${Boost_LIBRARIES})
  Message(${CMAKE_SOURCE_DIR}/tests/data)
  add_test (NAME ${filename} 
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests/data
    COMMAND ${filename})
ENDFOREACH(TEST_FILE ${TEST_FILES})

