if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
                message(STATUS "64-bit windows build detected, setting fftw dir to ${PROJECT_SOURCE_DIR}/fftw-3.2.2-dll64/")
		set( WIN_FFTW_LIBDIR "${PROJECT_SOURCE_DIR}/fftw-3.2.2-dll64/" )
	else( CMAKE_SIZEOF_VOID_P EQUAL 8 )
                message(STATUS "32-bit windows build detected, setting fftw dir to ${PROJECT_SOURCE_DIR}/fftw-3.2.2-dll32/")
		set( WIN_FFTW_LIBDIR "${PROJECT_SOURCE_DIR}/fftw-3.2.2-dll32/" )
	endif( CMAKE_SIZEOF_VOID_P EQUAL 8 )
endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

find_path(FFTW3_INCLUDE_DIRS fftw3.h HINTS $ENV{HOME}/include /usr/include "${WIN_FFTW_LIBDIR}")

IF(WIN32)
	find_library(FFTW3_LIBS libfftw3-3 "${WIN_FFTW_LIBDIR}" /usr/lib)
	find_library(FFTW3F_LIBS libfftw3f-3 "${WIN_FFTW_LIBDIR}" /usr/lib)
ELSEIF(UNIX)
	find_library(FFTW3_LIBS fftw3 HINTS $ENV{HOME}/lib /usr/lib)
	find_library(FFTW3F_LIBS fftw3f HINTS $ENV{HOME}/lib /usr/lib)
ENDIF(WIN32)

set(FFTW3_FOUND TRUE)
set(FFTW3F_FOUND TRUE)
 
if (NOT FFTW3_INCLUDE_DIRS)
  set(FFTW3_FOUND FALSE)
endif (NOT FFTW3_INCLUDE_DIRS)

if (NOT FFTW3_LIBS)
  set(FFTW3_FOUND FALSE)
endif (NOT FFTW3_LIBS)

if (NOT FFTW3F_LIBS)
  set(FFTW3F_FOUND FALSE)
endif (NOT FFTW3F_LIBS)
