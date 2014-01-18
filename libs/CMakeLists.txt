cmake_minimum_required(VERSION 2.8)

project(qstem)

#set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--export-all-symbols")
#set(CMAKE_SHARED_LINKER_FLAGS "-Wl")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC" )

#set(Boost_USE_STATIC_LIBS ON)

find_package(Boost COMPONENTS system filesystem REQUIRED)

FIND_PACKAGE(HDF5 "1.8.5" REQUIRED)

# Make boost, fftw3, and our headers available to all projects

FILE(GLOB LIBS_C_FILES "*.cpp")
FILE(GLOB LIBS_H_FILES "*.hpp")

# defines qstem_config (config readers)
add_subdirectory(config_IO)
# defines qstem_data (I/O for image data)
add_subdirectory(data_IO)
#defines qstem_structure
add_subdirectory(structure_IO)
# defines qstem_wavefunctions
add_subdirectory(wavefunctions)
# defines qstem_potentials
add_subdirectory(potentials)
# defines qstem_experiments
add_subdirectory(experiments)

include_directories(
  ${PROJECT_SOURCE_DIR}
  config_IO
  data_IO
  structure_IO
  wavefunctions
  potentials
  experiments
  ${Boost_INCLUDE_DIRS} 
  ${FFTW3_INCLUDE_DIRS} 
  ${HDF5_INCLUDE_DIR}  )

add_library(${PROJECT_NAME} SHARED ${LIBS_C_FILES})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} ${HDF5_LIBRARIES} ${FFTW_LIBRARIES})
