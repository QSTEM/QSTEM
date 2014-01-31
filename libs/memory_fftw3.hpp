/*
QSTEM - image simulation for TEM/STEM/CBED
    Copyright (C) 2000-2010  Christoph Koch
	Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MEMORY_H
#define MEMORY_H

// #include <stdlib.h>
#include <stdio.h>
// #include "floatdef.hpp"
#include "fftw3.h"
#include <boost/shared_array.hpp>


#include "stemtypes_fftw3.hpp"

/*---------------------------- float1D() -------------------------------*/
/*
	1D array allocator for type float
	make space for m[0...(n-1)]
	printf error message and exit if not successful
	
	this save checking for a NULL return etc every time 
	
*/
boost::shared_array<float_tt> float1D( int n, const char *message );

/*---------------------------- double1D() -------------------------------*/
/*
	1D array allocator for type double
	make space for m[0...(n-1)]
	printf error message and exit if not successful
	
	this save checking for a NULL return etc every time 
	
*/
boost::shared_array<double> double1D( int n, const char *message );


/*---------------------------- short2D() -------------------------------*/
/*
	2D array allocator for type short
	make space for m[0...(nx-1)][0..(ny-1)]
	
	message = char[] with error message
	
*/
short **short2D( int nx, int ny, const char *message );

/*---------------------------- long2D() -------------------------------*/
/*
	2D array allocator for type long
	make space for m[0...(nx-1)][0..(ny-1)]
	
	message = char[] with error message
	
*/
long **long2D( int nx, int ny, const char *message );
int **int2D( int nx, int ny, const char *message );

/*---------------------------- float2D() -------------------------------*/
/*
	2D array allocator for type float
	make space for m[0...(nx-1)][0..(ny-1)]

*/
float_tt **float2D( int nx, int ny, const char *message );
float **float32_2D( int nx, int ny, const char *message );
float ***float32_3D( int nx, int ny,int nz, const char *message );

/*---------------------------- double2D() -------------------------------*/
/*
	2D array allocator for type double
	make space for m[0...(nx-1)][0..(ny-1)]

*/
float_tt ***float3D( int nx, int ny,int nz, const char *message );
double **double2D( int nx, int ny, const char *message );

boost::shared_array<complex_tt> complex1D(int n, const char *message);
complex_tt  **complex2D(int nx, int ny, const char *message);
fftw_complex **complex2Dd(int nx, int ny, const char *message);
complex_tt  ***complex3D(int nx, int ny,int nz, const char *message);
fftw_complex  ***complex3Dd(int nx, int ny,int nz, const char *message);

void **any2D( int nx, int ny,int size, const char *message );
void ***any3D( int nx, int ny,int nz,int size, const char *message );

#endif










