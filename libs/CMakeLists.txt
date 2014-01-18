cmake_minimum_required(VERSION 2.8)

project(qstem)

#set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--export-all-symbols")
#set(CMAKE_SHARED_LINKER_FLAGS "-Wl")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC" )

#set(Boost_USE_STATIC_LIBS ON)

find_package(Boost COMPONENTS system filesystem REQUIRED)

FIND_PACKAGE(HDF5 "1.8.5" REQUIRED)

FILE(GLOB LIBS_FILES "*.?pp")

function (add_lib_dir dirname)
	#add_subdirectory(${dirname})
        #exploring this: not ideal to use glob because it defeats some of CMake's ability to detect changes
        #  a better way would be to manually add files in CMakeLists.txt in each folder (still adding them ultimately
        #  to LIBS_FILES
	FILE(GLOB TMP "${dirname}/*.?pp")
	set(LIBS_C_FILES ${LIBS_C_FILES} ${TMP} PARENT_SCOPE)
	SOURCE_GROUP(${dirname} FILES ${TMP} PARENT_SCOPE)
endfunction(add_lib_dir)

# defines qstem_config (config readers)
add_lib_dir(config_IO)
# defines qstem_data (I/O for image data)
add_lib_dir(data_IO)
#defines qstem_structure
add_lib_dir(structure_IO)
# defines qstem_wavefunctions
add_lib_dir(wavefunctions)
# defines qstem_potentials
add_lib_dir(potentials)
# defines qstem_experiments
add_lib_dir(experiments)

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

add_library(${PROJECT_NAME} SHARED ${LIBS_FILES})
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} ${HDF5_LIBRARIES} ${FFTW_LIBRARIES})
