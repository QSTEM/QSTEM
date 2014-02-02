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

FILE(GLOB QSTEM_LIB_HEADERS "${CMAKE_SOURCE_DIR}/libs/*.h")
FILE(GLOB STEM3_HEADERS "${CMAKE_SOURCE_DIR}/stem3/*.h" )
include_directories("${CMAKE_SOURCE_DIR}/libs" "${FFTW3_INCLUDE_DIRS}" ${Boost_INCLUDE_DIRS} "${CMAKE_SOURCE_DIR}/tests")
link_directories(${Boost_LIBRARY_DIRS})

#debugging - no QSTEM libraries at all.
#add_executable(test_libs test_main.cpp ${LIB_TEST_FILES})
#target_link_libraries(test_libs ${Boost_LIBRARIES})
add_subdirectory(libs)

#add_executable(test_stem3 test_main.cpp  ${STEM3_TEST_FILES} ${STEM3_TEST_HEADERS} ${STEM3_HEADERS} ${QSTEM_LIB_HEADERS})
#target_link_libraries(test_stem3 qstem_libs ${FFTW3_LIBS} ${FFTW3F_LIBS} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

#add_executable(test_gbmaker test_main.cpp  ${GBMAKER_TEST_FILES} ${GBMAKER_TEST_HEADERS} ${QSTEM_LIB_HEADERS})
#target_link_libraries(test_gbmaker qstem_libs ${FFTW3_LIBS} ${FFTW3F_LIBS} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

if (WIN32)
  #get_filename_component(UTF_RELEASE_NAME  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE} ABSOLUTE)
  #get_filename_component(UTF_DEBUG_NAME  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG} ABSOLUTE)
	
  get_filename_component(UTF_BASE_NAME 
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE} NAME_WE)
  get_filename_component(UTF_PATH 
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE} PATH)
  SET(UTF_RELEASE_NAME ${UTF_PATH}/${UTF_BASE_NAME}.dll)

  get_filename_component(UTF_BASE_NAME 
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG} NAME_WE)
  get_filename_component(UTF_PATH 
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG} PATH)
  SET(UTF_DEBUG_NAME ${UTF_PATH}/${UTF_BASE_NAME}.dll)

  add_custom_command(TARGET test_libs POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${UTF_RELEASE_NAME} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release)

  add_custom_command(TARGET test_libs POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    ${UTF_DEBUG_NAME} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug)
endif(WIN32)
