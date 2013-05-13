#ifndef MEMORY_H
#define MEMORY_H

// #include <stdlib.h>
#include <stdio.h>
// #include "floatdef.h"
#include "fftw3.h"

#ifndef float_tt
#define float_tt float
#endif


/*---------------------------- float1D() -------------------------------*/
/*
	1D array allocator for type float
	make space for m[0...(n-1)]
	printf error message and exit if not successful
	
	this save checking for a NULL return etc every time 
	
*/
//float_tt *float1D( int n, const char *message );

/*---------------------------- double1D() -------------------------------*/
/*
	1D array allocator for type double
	make space for m[0...(n-1)]
	printf error message and exit if not successful
	
	this save checking for a NULL return etc every time 
	
*/
//double* double1D( int n, const char *message );


/*---------------------------- short2D() -------------------------------*/
/*
	2D array allocator for type short
	make space for m[0...(nx-1)][0..(ny-1)]
	
	message = char[] with error message
	
*/
//short **short2D( int nx, int ny, const char *message );

/*---------------------------- long2D() -------------------------------*/
/*
	2D array allocator for type long
	make space for m[0...(nx-1)][0..(ny-1)]
	
	message = char[] with error message
	
*/
//long **long2D( int nx, int ny, const char *message );
//int **int2D( int nx, int ny, const char *message );

/*---------------------------- float2D() -------------------------------*/
/*
	2D array allocator for type float
	make space for m[0...(nx-1)][0..(ny-1)]

*/
//float_tt **float2D( int nx, int ny, const char *message );
//float **float32_2D( int nx, int ny, const char *message );
//float ***float32_3D( int nx, int ny,int nz, const char *message );

/*---------------------------- double2D() -------------------------------*/
/*
	2D array allocator for type double
	make space for m[0...(nx-1)][0..(ny-1)]

*/
//float_tt ***float3D( int nx, int ny,int nz, const char *message );
//double **double2D( int nx, int ny, const char *message );

//fftw_complex  **complex2D(int nx, int ny, const char *message);
//fftwf_complex **complex2Df(int nx, int ny, const char *message);  // single precision
//fftw_complex  ***complex3D(int nx, int ny,int nz, const char *message);
//fftwf_complex ***complex3Df(int nx, int ny,int nz, const char *message); // single precision

//void **any2D( int nx, int ny,int size, const char *message );
//void ***any3D( int nx, int ny,int nz,int size, const char *message );

#endif