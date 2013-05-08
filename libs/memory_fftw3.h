#ifndef MEMORY_H
#define MEMORY_H

// #include <stdlib.h>
#include <stdio.h>

#include "stemtypes_fftw3.h"

/************** Functions to get vectors/arrays out ****************/
boost::shared_ptr<float1D_type> float1D(int size, std::string message);
boost::shared_ptr<double1D_type> double1D(int size, std::string message);

boost::shared_ptr<float2D_type> float2D(int nx, int ny, std::string message);
boost::shared_ptr<double2D_type> double2D(int nx, int ny, std::string message);
boost::shared_ptr<int2D_type> int2D(int nx, int ny, std::string message);
boost::shared_ptr<complex2D_type> complex2D(int nx, int ny,int nz, const char *message);

boost::shared_ptr<float3D_type> float3D(int nx, int ny, std::string message);
boost::shared_ptr<complex3D_type> complex3D(int nx, int ny,int nz, const char *message);

void **any2D( int nx, int ny,int size, const char *message );
void ***any3D( int nx, int ny,int nz,int size, const char *message );

#endif

