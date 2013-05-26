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

#ifndef IMAGELIB_H
#define IMAGELIB_H

// #include "fftw.h"
// #include "floatdef.h"
#include "stemtypes_fftw3.h"


typedef struct imageStructType {
  int headerSize;  // first byte of image will be size of image header (in bytes)
                   // This is the size without the data and comment pointers!!!
  int paramSize;   // number of additional parameters
  int commentSize; // length of comment string
  int nx,ny;
  int complexFlag;
  int dataSize;    // size of one data element in bytes (e.g. complex double: 16)
  int version;     // The version flag will later help to find out how to 
                   // distinguish between images produced by different versions of stem
  double t;        // thickness
  double dx,dy;    // size of one pixel
  double *params;  // array for additional parameters
  char *comment;   // comment of prev. specified length
} imageStruct;


void getImageHeader(imageStruct *header,FILE * fp);
imageStruct *makeNewHeader(int nx,int ny);
imageStruct *makeNewHeaderCompact(int cFlag,int nx,int ny,double t,double dx,double dy,
				  int paramSize, double *params,char *comment); 
void setHeaderComment(imageStruct *header, char *comment);
  
imageStruct *readImage(void ***pix,int nx,int ny,char *fileName);
void writeImage(void **pix, imageStruct *header, char *fileName);
void writeRealImage(void **pix, imageStruct *header, char *fileName, int dataSize);
/*
void writeRealImage(fftw_real **pix, int nx, int ny, float_t dx, 
		   float_t dy, float_t t,char *fileName);
*/

void readRealImage(fftw_real **pix, int nx, int ny,real *dx, 
		   real *dy, real *t, char *fileName);

// old image I/O functions:
// void readRealImage_old(fftw_real **pix, int nx, int ny,float_t *t, char *fileName);
// void readImage_old(fftw_complex **pix, int nx, int ny,float_t *t, char *fileName);
// void writeRealImage_old(fftw_real **pix, int nx, int ny, float_t t,char *fileName);
// void writeImage_old(fftw_complex **pix, int nx, int ny, float_t t,char *fileName);


/**************************************************************
 * Here is how to use the new image writing routines
 *
 * static imageStruct *header = NULL;
 *
 * if (header == NULL) header = makeNewHeaderCompact(cFlag,Nx,Ny,t,dx,dy,0,NULL,comment);
 * writeImage(cimage,header,filename);
 * or : writeRealImage(rimage,header,filename,sizeof(float));
 *
 **************************************************************
 * Reading an image works like this:
 * 
 * imageStruct *header;
 * header = readImage((void ***)(&pix),nx,ny,fileName);
 *
 * [This function will read an image.  It reuses the same header 
 * struct over and over.  Therefore, values must be copied from 
 * the header members before calling this function again.
 *
 * The image pointer may also be NULL, in which case memory will be
 * allocated for it, and its size will be returned in the header struct
 * members nx, and ny.]
 **************************************************************/

#endif
