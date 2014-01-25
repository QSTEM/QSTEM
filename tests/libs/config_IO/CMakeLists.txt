cmake_minimum_required(VERSION 2.8.3)

PROJECT(Test_configIO)

FILE(GLOB TEST_FILES "*.cpp")
FOREACH(TEST_FILE ${TEST_FILES})
  GET_FILENAME_COMPONENT(filename ${TEST_FILE} NAME_WE)
  add_executable(${filename} ${TEST_FILE})
  target_link_libraries(${filename} qstem ${FFTW3_LIBS} ${FFTW3F_LIBS} ${Boost_LIBRARIES})
  add_test (NAME ${filename} 
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/tests/data"
    COMMAND ${filename})
ENDFOREACH(TEST_FILE ${TEST_FILES})
