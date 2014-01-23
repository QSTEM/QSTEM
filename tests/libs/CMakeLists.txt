cmake_minimum_required(VERSION 2.8)

project(tests)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test_bin)

add_definitions (-DBOOST_TEST_DYN_LINK)
add_definitions(-DBOOST_ALL_NO_LIB)  # tell the compiler to undefine this boost macro

set(Boost_USE_STATIC_LIBS        OFF)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     OFF)
set(BOOST_ALL_DYN_LINK           ON)   # force dynamic linking for all libraries

find_package(Boost 1.35.0 COMPONENTS unit_test_framework REQUIRED)

include_directories("${CMAKE_SOURCE_DIR}/libs" "${FFTW3_INCLUDE_DIRS}" ${Boost_INCLUDE_DIRS})

#debugging - no QSTEM libraries at all.
#add_executable(test_libs test_main.cpp ${LIB_TEST_FILES})
#target_link_libraries(test_libs ${Boost_LIBRARIES})

FILE(GLOB CONFIG_TESTS "config_IO/*.cpp")
add_executable(test_config_reader ${CONFIG_TESTS})
target_link_libraries(test_config_reader qstem ${FFTW3_LIBS} ${FFTW3F_LIBS} ${Boost_LIBRARIES})
add_test (NAME TestConfigReader WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/tests/data"
  COMMAND test_config_reader)


add_executable(test_wavefunctions test_wavefunctions.cpp)
target_link_libraries(test_wavefunctions qstem ${FFTW3_LIBS} ${FFTW3F_LIBS} ${Boost_LIBRARIES})
add_test (NAME TestWaveFunctions WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/tests/data"
  COMMAND test_wavefunctions)


