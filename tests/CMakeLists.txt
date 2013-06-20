cmake_minimum_required(VERSION 2.8)

project(tests)

set (Boost_USE_STATIC_LIBS ON)
add_definitions( -DBOOST_ALL_NO_LIB )

FILE(GLOB QSTEM_LIB_HEADERS "${CMAKE_SOURCE_DIR}/libs/*.h")
FILE(GLOB STEM3_HEADERS "${CMAKE_SOURCE_DIR}/stem3/*.h")
include_directories("${CMAKE_SOURCE_DIR}/libs" "${FFTW3_INCLUDE_DIRS}")
link_directories(${Boost_LIBRARY_DIRS})

FILE(GLOB STEM3_TEST_HEADERS "${CMAKE_SOURCE_DIR}/stem3/tests/*.h")
FILE(GLOB STEM3_TEST_FILES "${CMAKE_SOURCE_DIR}/stem3/tests/*.cpp")
FILE(GLOB LIB_TEST_HEADERS "${CMAKE_SOURCE_DIR}/libs/tests/*.h")
FILE(GLOB LIB_TEST_FILES "${CMAKE_SOURCE_DIR}/libs/tests/*.cpp")
FILE(GLOB GBMAKER_TEST_HEADERS "${CMAKE_SOURCE_DIR}/gbmaker/tests/*.h")
FILE(GLOB GBMAKER_TEST_FILES "${CMAKE_SOURCE_DIR}/gbmaker/tests/*.cpp")

add_definitions (-DBOOST_TEST_DYN_LINK)

#add_executable(test_libs ${LIB_TEST_FILES} ${LIB_TEST_HEADERS} ${QSTEM_LIB_HEADERS})
add_executable(test_libs ${LIB_TEST_FILES})
#target_link_libraries(test_libs qstem_libs ${FFTW3_LIBS} ${FFTW3F_LIBS} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
target_link_libraries(test_libs ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

#add_executable(test_stem3 ${STEM3_TEST_FILES} ${STEM3_TEST_HEADERS} ${STEM3_HEADERS} ${QSTEM_LIB_HEADERS})
#target_link_libraries(test_stem3 qstem_libs ${FFTW3_LIBS} ${FFTW3F_LIBS} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

#add_executable(test_gbmaker ${GBMAKER_TEST_FILES} ${GBMAKER_TEST_HEADERS} ${QSTEM_LIB_HEADERS})
#target_link_libraries(test_gbmaker qstem_libs ${FFTW3_LIBS} ${FFTW3F_LIBS} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
