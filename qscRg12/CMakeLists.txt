cmake_minimum_required(VERSION 2.8)

project(qscRg12)

FILE(GLOB QSCRG12 "${CMAKE_SOURCE_DIR}/qscRg12/*.c")

add_executable(qscRg12 ${QSCRG12})
target_link_libraries(qscRg12 ${M_LIB}) 
