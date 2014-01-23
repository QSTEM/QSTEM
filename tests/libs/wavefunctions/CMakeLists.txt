cmake_minimum_required(VERSION 2.8.3)

PROJECT(Test_Wavefunctions)

FILE(GLOB TEST_FILES "*.cpp")
add_executable(Test_Wavefunctions ${TEST_FILES})
target_link_libraries(${PROJECT_NAME} qstem ${FFTW3_LIBS} ${FFTW3F_LIBS} ${Boost_LIBRARIES})
add_test (NAME Test_Wave WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/tests/data"
	  COMMAND Test_Wavefunctions)