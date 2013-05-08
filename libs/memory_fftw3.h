#ifndef MEMORY_H
#define MEMORY_H

// #include <stdlib.h>
#include <stdio.h>

#include "stemtypes_fftw3.h"

/************** Functions to get vectors/arrays out ****************/
float1DArray float1D(int size, std::string message);
double1DArray double1D(int size, std::string message);
int1DArray int1D(int size, std::string message);

float2DArray float2D(int nx, int ny, std::string message);
double2DArray double2D(int nx, int ny, std::string message);
int2DArray int2D(int nx, int ny, std::string message);
complex2DArray complex2D(int nx, int ny, std::string message);

float3DArray float3D(int nx, int ny, int nz, std::string message);
complex3DArray complex3D(int nx, int ny,int nz, std::string message);

void **any2D( int nx, int ny,int size, const char *message );
void ***any3D( int nx, int ny,int nz,int size, const char *message );

#endif

